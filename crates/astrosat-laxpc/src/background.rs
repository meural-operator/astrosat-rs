//! Stage 6: Background modeling (2D/3D B-spline + 5-harmonic diurnal variation),
//! `gaina8.dat` / `lxhv.dat` gain tracking, background spectral & light-curve shifting,
//! and automatic RMF recommendation.
//!
//! Faithfully implements the physics of `laxpcl1.f` (`BSPLIN`, `BSPEV2`, `BSPEVN`, diurnal harmonics)
//! and `backshiftv3.f` (`gainsh`, HV-jump-aware interpolation, secular polynomial scaling, RMF selection).

use crate::error::{LaxpcError, Result};
use crate::types::{AnodeSelection, LaxpcUnit};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Default diurnal background oscillation period (~84496 s, `readme.laxpcl1:53`; stored in `backfit.*` as `84494.9 s`).
pub const DEFAULT_DIURNAL_PERIOD_SEC: f64 = 84_496.0;

/// Channel-to-energy table loaded from `lxp{10,20,30}.channel`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelMap {
    /// Lower energy edge (keV) for each channel `0..nch`.
    pub e_min_kev: Vec<f64>,
    /// Upper energy edge (keV) for each channel `0..nch`.
    pub e_max_kev: Vec<f64>,
}

impl ChannelMap {
    /// Parse an ASCII `lxp{10,20,30}.channel` file (two columns `e_min e_max` in keV).
    pub fn from_file(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|e| LaxpcError::ChannelMap {
            path: path.to_path_buf(),
            message: e.to_string(),
        })?;
        let reader = BufReader::new(file);
        let mut e_min_kev = Vec::with_capacity(512);
        let mut e_max_kev = Vec::with_capacity(512);

        for (line_no, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(LaxpcError::ChannelMap {
                    path: path.to_path_buf(),
                    message: format!("line {}: expected 2 columns, got {}", line_no + 1, parts.len()),
                });
            }
            let lo: f64 = parts[0].parse().map_err(|_| LaxpcError::ChannelMap {
                path: path.to_path_buf(),
                message: format!("line {}: invalid e_min '{}'", line_no + 1, parts[0]),
            })?;
            let hi: f64 = parts[1].parse().map_err(|_| LaxpcError::ChannelMap {
                path: path.to_path_buf(),
                message: format!("line {}: invalid e_max '{}'", line_no + 1, parts[1]),
            })?;
            e_min_kev.push(lo);
            e_max_kev.push(hi);
        }

        if e_min_kev.is_empty() {
            return Err(LaxpcError::ChannelMap {
                path: path.to_path_buf(),
                message: "empty channel map".to_string(),
            });
        }

        Ok(Self {
            e_min_kev,
            e_max_kev,
        })
    }

    /// Return the nominal or post-anomaly shifted `.channel` filename for `unit` at mission day
    /// `ndays`, matching `backshiftv3.f:334-346`.
    pub fn filename_for_day(unit: LaxpcUnit, ndays: f64) -> &'static str {
        match unit {
            LaxpcUnit::Lxp10 => {
                if ndays > 1112.0 {
                    "lxp10cshm212v2.2.channel"
                } else if ndays > 1035.0 {
                    "lxp10cshm197v2.2.channel"
                } else if ndays > 990.0 {
                    "lxp10cshm175v2.2.channel"
                } else if ndays > 917.0 {
                    "lxp10cshm141v2.2.channel"
                } else {
                    "lxp10.channel"
                }
            }
            LaxpcUnit::Lxp20 => "lxp20.channel",
            LaxpcUnit::Lxp30 => "lxp30.channel",
        }
    }

    /// Load the epoch-appropriate `.channel` map from `caldb_dir` for `unit` at mission day `ndays`
    /// (`backshiftv3.f:334-346`), falling back to `lxp{10,20,30}.channel` if the shifted file is absent.
    pub fn from_caldb_for_day(caldb_dir: &Path, unit: LaxpcUnit, ndays: f64) -> Option<Self> {
        let preferred = caldb_dir.join(Self::filename_for_day(unit, ndays));
        if preferred.is_file() {
            if let Ok(cmap) = Self::from_file(&preferred) {
                return Some(cmap);
            }
        }
        let fallback = caldb_dir.join(format!("{}0.channel", unit.prefix()));
        if fallback != preferred && fallback.is_file() {
            return Self::from_file(&fallback).ok();
        }
        None
    }

    /// Number of channels in this map (`512` for Lxp10/30, `256` for Lxp20).
    #[inline]
    pub fn len(&self) -> usize {
        self.e_min_kev.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.e_min_kev.is_empty()
    }

    /// Convert an energy range `(e_low_kev, e_high_kev)` into inclusive `(chan_min, chan_max)`
    /// in native channel indices (`0..len()`).
    pub fn energy_to_native_channels(&self, e_low_kev: f64, e_high_kev: f64) -> (usize, usize) {
        let mut lo = self.len().saturating_sub(1);
        let mut hi = 0usize;
        for i in 0..self.len() {
            if self.e_max_kev[i] >= e_low_kev && i < lo {
                lo = i;
            }
            if self.e_min_kev[i] <= e_high_kev {
                hi = i;
            }
        }
        if lo > hi {
            (0, self.len().saturating_sub(1))
        } else {
            (lo, hi)
        }
    }

    /// Convert an energy range `(e_low_kev, e_high_kev)` into inclusive 1024-scale raw PHA channel bounds
    /// `(pha_min, pha_max)` suitable for `EventFilter`.
    pub fn energy_to_pha_1024(&self, e_low_kev: f64, e_high_kev: f64) -> (u16, u16) {
        let (lo, hi) = self.energy_to_native_channels(e_low_kev, e_high_kev);
        let scale = (1024 / self.len().max(1)).max(1);
        let pha_min = (lo * scale).min(1023) as u16;
        let pha_max = ((hi + 1) * scale - 1).min(1023) as u16;
        (pha_min, pha_max)
    }
}

/// One row of `gaina8.dat` (days since launch + 30 keV / 60 keV peak & FWHM for all 3 LAXPC units).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GainEntry {
    /// Days since AstroSat launch (`ndays`).
    pub day: f64,
    /// 30 keV calibration peak channel in anode A8 for `[Lxp10, Lxp20, Lxp30]`.
    pub peak_30kev: [f64; 3],
    /// 30 keV FWHM resolution for `[Lxp10, Lxp20, Lxp30]`.
    pub fwhm_30kev: [f64; 3],
    /// 60 keV calibration peak channel in anode A8 for `[Lxp10, Lxp20, Lxp30]`.
    pub peak_60kev: [f64; 3],
    /// 60 keV FWHM resolution for `[Lxp10, Lxp20, Lxp30]`.
    pub fwhm_60kev: [f64; 3],
}

/// Parsed `gaina8.dat` gain history table.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GainTable {
    pub entries: Vec<GainEntry>,
}

impl GainTable {
    /// Parse a `gaina8.dat` file (skipping `#` comments and blank lines).
    pub fn from_file(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(BufReader::new(file))
    }

    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        let mut entries = Vec::new();
        for line_res in reader.lines() {
            let line = line_res?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            // Fortran list-directed rows after LAXPC30 shutoff end with `/` and have 9 numeric columns
            let line_no_slash = trimmed.split('/').next().unwrap_or(trimmed);
            let vals: Vec<f64> = line_no_slash
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if vals.len() < 9 {
                continue;
            }
            let p30_3 = vals.get(9).copied().unwrap_or(0.0);
            let f30_3 = vals.get(10).copied().unwrap_or(0.0);
            let p60_3 = vals.get(11).copied().unwrap_or(0.0);
            let f60_3 = vals.get(12).copied().unwrap_or(0.0);
            entries.push(GainEntry {
                day: vals[0],
                peak_30kev: [vals[1], vals[5], p30_3],
                fwhm_30kev: [vals[2], vals[6], f30_3],
                peak_60kev: [vals[3], vals[7], p60_3],
                fwhm_60kev: [vals[4], vals[8], f60_3],
            });
        }
        if entries.is_empty() {
            return Err(LaxpcError::Background(
                "gaina8.dat contained no valid gain entries".to_string(),
            ));
        }
        entries.sort_by(|a, b| a.day.total_cmp(&b.day));
        Ok(Self { entries })
    }

    /// Maximum mission day tabulated in `gaina8.dat`.
    pub fn max_day(&self) -> f64 {
        self.entries.last().map(|e| e.day).unwrap_or(0.0)
    }

    /// Interpolate 30 keV peak channel at `ndays` for `unit`, respecting HV jump epochs in `hv_table`
    /// (`backshiftv3.f:248-307`).
    pub fn peak_30kev_at(&self, unit: LaxpcUnit, ndays: f64, hv_table: Option<&HvJumpTable>) -> Result<f64> {
        let u = unit.index();
        let days: Vec<f64> = self.entries.iter().map(|e| e.day).collect();
        let peaks: Vec<f64> = self.entries.iter().map(|e| e.peak_30kev[u]).collect();

        if unit == LaxpcUnit::Lxp30 && ndays > 892.0 {
            return Err(LaxpcError::Background(format!(
                "LAXPC30 was switched off after mission day 892 (requested day {ndays:.1})"
            )));
        }

        // Standard 2-point linear interpolation (`DIVDIF` with `NUSE=2`)
        let mut chs = crate::bb::interp_linear(&days, &peaks, ndays);

        // One-sided interpolation near HV jump epochs (`backshiftv3.f:258-292`)
        if let Some(hv) = hv_table {
            if let Some(tgain) = hv.nearest_jump_within(unit, ndays, 1.0) {
                let mut it = days.partition_point(|&d| d < ndays);
                if it >= days.len() {
                    it = days.len() - 1;
                } else if days[it] > ndays && it > 0 {
                    it -= 1;
                }
                if tgain < ndays {
                    if it + 1 < days.len() {
                        let dt = days[it + 1] - days[it];
                        chs = if dt.abs() > 0.0 {
                            peaks[it] + (ndays - days[it]) * (peaks[it + 1] - peaks[it]) / dt
                        } else {
                            peaks[it]
                        };
                    } else {
                        chs = peaks[it];
                    }
                } else if tgain > ndays && it > 0 {
                    let dt = days[it] - days[it - 1];
                    chs = if dt.abs() > 0.0 {
                        peaks[it] + (ndays - days[it]) * (peaks[it] - peaks[it - 1]) / dt
                    } else {
                        peaks[it]
                    };
                }
            }
        }

        // Post-day-883 LAXPC30 exponential gain extrapolation (`backshiftv3.f:294`)
        if unit == LaxpcUnit::Lxp30 && ndays > 883.0 {
            chs = 450.0 * ((ndays - 871.0) / 30.0).exp();
        }

        Ok(chs)
    }
}

/// Parsed `lxhv.dat` HV adjustment epoch table.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct HvJumpTable {
    /// Per-unit sorted list of mission days when HV was adjusted (`[Lxp10, Lxp20, Lxp30]`).
    pub jumps_by_unit: [Vec<f64>; 3],
}

impl HvJumpTable {
    /// Parse an `lxhv.dat` file (two columns: `day unit_id`).
    pub fn from_file(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(BufReader::new(file))
    }

    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        let mut table = Self::default();
        for line_res in reader.lines() {
            let line = line_res?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() < 2 {
                continue;
            }
            if let (Ok(day), Ok(uid)) = (parts[0].parse::<f64>(), parts[1].parse::<u8>()) {
                if let Ok(unit) = LaxpcUnit::from_id(uid) {
                    table.jumps_by_unit[unit.index()].push(day);
                }
            }
        }
        for v in &mut table.jumps_by_unit {
            v.sort_by(|a, b| a.total_cmp(b));
        }
        Ok(table)
    }

    /// Return the jump epoch if `ndays` lies within `tol_days` of any HV adjustment for `unit`.
    pub fn nearest_jump_within(&self, unit: LaxpcUnit, ndays: f64, tol_days: f64) -> Option<f64> {
        self.jumps_by_unit[unit.index()]
            .iter()
            .copied()
            .find(|&th| (th - ndays).abs() <= tol_days)
    }

    /// Return all HV adjustment epochs strictly between `day_a` and `day_b` for `unit`.
    pub fn jumps_between(&self, unit: LaxpcUnit, day_a: f64, day_b: f64) -> Vec<f64> {
        let lo = day_a.min(day_b);
        let hi = day_a.max(day_b);
        self.jumps_by_unit[unit.index()]
            .iter()
            .copied()
            .filter(|&th| th > lo && th < hi)
            .collect()
    }
}

/// One entry in the TIFR `backshiftv3.inp` blank-sky background catalog.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BackgroundCatalogEntry {
    /// Short epoch tag (e.g. `"jan16"`, `"aug17"`, `"apr18"`).
    pub short_name: String,
    /// Mission day since AstroSat launch (`ndy`).
    pub day: f64,
    /// Measured 30 keV calibration peak channel in anode A8 for `[Lxp10, Lxp20, Lxp30]` (`999.0` = unit off).
    pub peak_30kev: [f64; 3],
    /// Nominal LAXPC10 background `.pha` filename (e.g. `"backlxp10aug17.pha"`).
    pub filename: String,
}

/// Official TIFR LaxpcSoft v3.4.6 blank-sky background catalog (`backshiftv3.inp`),
/// mapping background epochs/filenames to their measured 30 keV A8 calibration peaks (`backshiftv3.f:108-165,295`).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BackgroundCatalog {
    pub entries: Vec<BackgroundCatalogEntry>,
}

/// Built-in fallback catalog matching `data/ref_src/laxpc/backshiftv3.inp` (used when `backshiftv3.inp`
/// is not present in `caldb_dir`, e.g. `data/caldb/laxpc/`).
const BUILTIN_BACKSHIFT_CATALOG: &[(&str, f64, [f64; 3])] = &[
    ("oct15", 24.0, [281.6, 221.3, 258.7]),
    ("nov15", 56.0, [255.7, 221.4, 251.2]),
    ("jan16", 100.0, [237.3, 223.6, 245.0]),
    ("mar16", 178.0, [247.6, 217.0, 232.7]),
    ("aug16", 323.0, [221.6, 211.9, 250.7]),
    ("sep16", 354.0, [228.4, 211.1, 265.5]),
    ("oct16", 384.0, [234.5, 209.5, 263.8]),
    ("dec16", 433.0, [243.6, 219.4, 233.1]),
    ("jan17", 483.0, [238.3, 217.3, 245.0]),
    ("feb17", 504.0, [239.2, 216.7, 238.6]),
    ("apr17", 563.0, [232.9, 228.6, 239.3]),
    ("may17", 610.0, [240.9, 226.9, 290.3]),
    ("aug17", 694.0, [234.4, 225.0, 213.0]),
    ("nov17", 785.0, [224.0, 224.7, 249.3]),
    ("jan18", 831.0, [230.5, 225.1, 247.5]),
    ("feb18", 870.0, [237.1, 222.5, 438.0]),
    ("apr18", 917.0, [157.8, 222.9, 999.0]),
    ("jun18", 982.0, [106.1, 222.4, 999.0]),
    ("jul18", 1012.0, [66.8, 222.1, 999.0]),
    ("jum18", 1021.0, [67.1, 224.1, 999.0]),
    ("aug18", 1039.0, [42.0, 222.5, 999.0]),
    ("auh18", 1045.0, [42.3, 223.1, 999.0]),
    ("aui18", 1049.0, [42.6, 223.7, 999.0]),
    ("auj18", 1052.0, [42.6, 223.5, 999.0]),
    ("oct18", 1120.0, [27.7, 223.3, 999.0]),
    ("nov18", 1151.0, [27.8, 221.2, 999.0]),
    ("feb19", 1231.0, [29.7, 221.2, 999.0]),
    ("mar19", 1267.0, [30.2, 223.2, 999.0]),
    ("apr19", 1300.0, [23.0, 220.6, 999.0]),
    ("may19", 1334.0, [23.0, 221.0, 999.0]),
    ("jun19", 1352.0, [23.0, 221.0, 999.0]),
    ("juo19", 1363.0, [23.0, 222.0, 999.0]),
    ("aug19", 1413.0, [23.0, 221.2, 999.0]),
    ("sep19", 1444.0, [23.0, 220.2, 999.0]),
    ("nov19", 1516.0, [23.9, 220.8, 999.0]),
    ("dec19", 1542.0, [24.1, 221.1, 999.0]),
    ("jan20", 1573.0, [24.5, 221.9, 999.0]),
    ("feb20", 1600.0, [24.9, 219.7, 999.0]),
    ("mar20", 1628.0, [25.5, 220.5, 999.0]),
    ("apr20", 1664.0, [26.2, 219.6, 999.0]),
    ("may20", 1696.0, [26.5, 218.6, 999.0]),
    ("jun20", 1721.0, [27.0, 219.0, 999.0]),
    ("aug20", 1771.0, [28.2, 218.0, 999.0]),
    ("auh20", 1792.0, [28.5, 217.7, 999.0]),
    ("sep20", 1821.0, [29.0, 217.0, 999.0]),
    ("nov20", 1861.0, [30.0, 217.9, 999.0]),
    ("dec20", 1907.0, [31.0, 217.4, 999.0]),
    ("jan21", 1944.0, [31.4, 217.4, 999.0]),
    ("feb21", 1962.0, [31.7, 219.0, 999.0]),
    ("mar21", 1997.0, [32.4, 217.6, 999.0]),
    ("apr21", 2034.0, [33.2, 217.6, 999.0]),
    ("may21", 2061.0, [34.2, 218.0, 999.0]),
    ("jun21", 2098.0, [34.7, 216.0, 999.0]),
    ("jul21", 2109.0, [34.9, 218.0, 999.0]),
    ("aug21", 2134.0, [35.4, 216.0, 999.0]),
    ("auh21", 2158.0, [35.7, 215.0, 999.0]),
    ("sep21", 2167.0, [36.1, 215.0, 999.0]),
    ("seq21", 2192.0, [36.6, 215.0, 999.0]),
    ("oct21", 2218.0, [37.0, 214.0, 999.0]),
    ("nov21", 2226.0, [37.7, 215.4, 999.0]),
    ("dec21", 2257.0, [38.0, 216.4, 999.0]),
    ("jan22", 2301.0, [39.0, 215.1, 999.0]),
    ("feb22", 2319.0, [40.0, 214.3, 999.0]),
    ("mar22", 2346.0, [40.0, 214.9, 999.0]),
    ("apr22", 2376.0, [41.0, 214.7, 999.0]),
    ("may22", 2409.0, [41.0, 214.5, 999.0]),
    ("jul22", 2479.0, [43.0, 214.0, 999.0]),
    ("aug22", 2497.0, [43.0, 213.9, 999.0]),
    ("auh22", 2526.0, [27.0, 213.7, 999.0]),
    ("sep22", 2532.0, [27.0, 213.6, 999.0]),
    ("seq22", 2537.0, [27.0, 213.7, 999.0]),
    ("oct22", 2557.0, [27.0, 213.5, 999.0]),
    ("ocu22", 2576.0, [28.0, 213.4, 999.0]),
    ("nov22", 2588.0, [28.0, 213.3, 999.0]),
    ("dec22", 2636.0, [28.0, 213.0, 999.0]),
    ("jan23", 2667.0, [29.0, 212.8, 999.0]),
    ("feb23", 2699.0, [30.0, 212.6, 999.0]),
    ("mar23", 2717.0, [30.0, 212.5, 999.0]),
    ("apr23", 2754.0, [31.0, 212.2, 999.0]),
    ("may23", 2784.0, [31.0, 212.0, 999.0]),
    ("jun23", 2809.0, [31.0, 211.9, 999.0]),
    ("juo23", 2827.0, [32.0, 211.8, 999.0]),
    ("jul23", 2847.0, [32.0, 211.7, 999.0]),
    ("jum23", 2859.0, [32.0, 211.6, 999.0]),
    ("aug23", 2869.0, [32.0, 211.5, 999.0]),
    ("auh23", 2875.0, [33.0, 211.4, 999.0]),
    ("aui23", 2892.0, [33.0, 211.3, 999.0]),
    ("sep23", 2909.0, [33.0, 211.2, 999.0]),
    ("oct23", 2945.0, [34.0, 211.0, 999.0]),
    ("nov23", 2964.0, [35.0, 210.9, 999.0]),
    ("jan24", 3025.0, [37.0, 227.9, 999.0]),
    ("mar24", 3083.0, [38.0, 227.6, 999.0]),
    ("apr24", 3114.0, [38.0, 227.4, 999.0]),
    ("may24", 3146.0, [39.0, 227.2, 999.0]),
    ("maz24", 3162.0, [40.0, 227.1, 999.0]),
    ("jul24", 3226.0, [42.0, 226.9, 999.0]),
    ("aug24", 3253.0, [43.0, 226.6, 999.0]),
    ("auh24", 3256.0, [43.0, 226.6, 999.0]),
    ("feb25", 3431.0, [48.0, 225.6, 999.0]),
    ("mar25", 3459.0, [49.0, 225.4, 999.0]),
    ("oct25", 3661.0, [54.0, 224.3, 999.0]),
    ("apr26", 3854.0, [59.0, 223.2, 999.0]),
    ("jul26", 3943.0, [61.0, 222.7, 999.0]),
];

impl BackgroundCatalog {
    /// Return the built-in TIFR LaxpcSoft v3.4.6 background catalog.
    pub fn builtin() -> Self {
        let entries = BUILTIN_BACKSHIFT_CATALOG
            .iter()
            .map(|&(name, day, peak_30kev)| BackgroundCatalogEntry {
                short_name: name.to_string(),
                day,
                peak_30kev,
                filename: format!("backlxp10{name}.pha"),
            })
            .collect();
        Self { entries }
    }

    /// Parse a `backshiftv3.inp` file (`backshiftv3.f:108-115`).
    pub fn from_file(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(BufReader::new(file))
    }

    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        let mut entries = Vec::new();
        for (idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            // First 9 lines of `backshiftv3.inp` are filenames + secular polynomial coefficients
            if idx < 9 {
                continue;
            }
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            // Expected: short_name obs_num dd mm yy nday g1 g2 g3 filename (10 tokens)
            if parts.len() < 10 {
                continue;
            }
            let (Ok(day), Ok(g1), Ok(g2), Ok(g3)) = (
                parts[5].parse::<f64>(),
                parts[6].parse::<f64>(),
                parts[7].parse::<f64>(),
                parts[8].parse::<f64>(),
            ) else {
                continue;
            };
            entries.push(BackgroundCatalogEntry {
                short_name: parts[0].to_ascii_lowercase(),
                day,
                peak_30kev: [g1, g2, g3],
                filename: parts[9].to_string(),
            });
        }
        if entries.is_empty() {
            return Err(LaxpcError::Background(
                "backshiftv3.inp contained no valid background catalog entries".to_string(),
            ));
        }
        Ok(Self { entries })
    }

    /// Load `backshiftv3.inp` from the first directory in `search_dirs` that contains it,
    /// falling back to [`BackgroundCatalog::builtin()`] when none is found on disk.
    pub fn load_or_builtin(search_dirs: &[Option<&Path>]) -> Self {
        for dir in search_dirs.iter().flatten() {
            let candidate = dir.join("backshiftv3.inp");
            if candidate.is_file() {
                if let Ok(cat) = Self::from_file(&candidate) {
                    return cat;
                }
            }
        }
        Self::builtin()
    }

    /// Look up the measured 30 keV calibration peak channel (`chb`) for a background spectrum
    /// by filename tag (e.g. `backlxp30aug17.pha` -> `"aug17"`) or by mission day (`|ndayb - ndy| < 2`,
    /// matching `backshiftv3.f:160-165,295`). Ignores `>= 999.0` off-unit sentinels.
    pub fn lookup_peak_30kev(
        &self,
        unit: LaxpcUnit,
        bkg_pha_path: Option<&Path>,
        bkg_day: f64,
    ) -> Option<f64> {
        let u = unit.index();

        // 1. Exact tag match against filename stem (e.g. `backlxp30aug17.pha`, `backlxp10apr18SE.pha`)
        if let Some(stem) = bkg_pha_path
            .and_then(|p| p.file_stem())
            .and_then(|s| s.to_str())
        {
            let lower = stem.to_ascii_lowercase();
            let stripped = lower.strip_suffix("se").unwrap_or(&lower);
            for entry in &self.entries {
                if stripped.ends_with(&entry.short_name) {
                    let peak = entry.peak_30kev[u];
                    if peak.is_finite() && peak > 0.0 && peak < 990.0 {
                        return Some(peak);
                    }
                }
            }
        }

        // 2. Mission-day match `|bkg_day - entry.day| < 2.0` (`backshiftv3.f:161`)
        //    (pick the last matching entry just like Fortran's loop `do i=1,nb`)
        let mut matched: Option<f64> = None;
        for entry in &self.entries {
            if (bkg_day - entry.day).abs() < 2.0 {
                let peak = entry.peak_30kev[u];
                if peak.is_finite() && peak > 0.0 && peak < 990.0 {
                    matched = Some(peak);
                }
            }
        }
        matched
    }
}

/// Result of computing gain shift between a background observation and a source observation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GainShiftSolution {
    pub unit: LaxpcUnit,
    pub source_day: f64,
    pub background_day: f64,
    /// Estimated 30 keV peak channel at source epoch.
    pub source_peak_ch: f64,
    /// 30 keV peak channel at background epoch.
    pub background_peak_ch: f64,
    /// Channel shift `dch = source_peak_ch - background_peak_ch` at 30 keV (1024-channel scale).
    pub delta_ch_30kev: f64,
    /// Relative gain scale factor `source_peak_ch / background_peak_ch`.
    pub gain_ratio: f64,
    /// Secular background count-rate trend scaling factor `poly(t_src)/poly(t_bkg)`.
    pub secular_scale: f64,
    /// HV jumps detected between `background_day` and `source_day` (non-empty = user warning).
    pub crossed_hv_jumps: Vec<f64>,
    /// Recommended XSPEC gain offset (keV) for LAXPC20 (`0.08919 - 0.0003001 * ndays`).
    pub lxp20_gain_offset_kev: Option<f64>,
    /// Recommended RMF filename matching TIFR CALDB conventions.
    pub recommended_rmf: String,
}

/// Default secular background polynomial coefficients from `backshiftv3.inp` (`maxday = 2321`).
pub const DEFAULT_SECULAR_MAXDAY_LXP20: f64 = 2321.0;
pub const DEFAULT_SECULAR_MAXDAY_OTHER: f64 = 875.0;
pub const DEFAULT_SECULAR_POLY_LXP10: &[f64] = &[269.273, -133.573, 1105.366, -2152.542, 1245.840];
pub const DEFAULT_SECULAR_POLY_LXP20: &[f64] = &[
    194.648, -71.6434, 524.2829, -1589.127, 2567.190, -2324.871, 1183.798, -316.5097, 34.55194,
];
pub const DEFAULT_SECULAR_POLY_LXP30: &[f64] = &[190.480, 24.874, 42.755, -210.185, 163.948];

fn eval_poly(coeffs: &[f64], x: f64) -> f64 {
    let mut acc = 0.0;
    for &c in coeffs.iter().rev() {
        acc = acc * x + c;
    }
    acc
}

/// Compute the secular background count-rate scaling factor `ss(source_day) / ss(background_day)`
/// matching `backshiftv3.f:325-333`.
pub fn secular_background_scale(unit: LaxpcUnit, source_day: f64, background_day: f64) -> f64 {
    let (maxday, coeffs) = match unit {
        LaxpcUnit::Lxp10 => (DEFAULT_SECULAR_MAXDAY_OTHER, DEFAULT_SECULAR_POLY_LXP10),
        LaxpcUnit::Lxp20 => (DEFAULT_SECULAR_MAXDAY_LXP20, DEFAULT_SECULAR_POLY_LXP20),
        LaxpcUnit::Lxp30 => (DEFAULT_SECULAR_MAXDAY_OTHER, DEFAULT_SECULAR_POLY_LXP30),
    };
    let t1 = source_day.min(maxday) / 1000.0;
    let t2 = background_day.min(maxday) / 1000.0;
    let s1 = eval_poly(coeffs, t1);
    let s2 = eval_poly(coeffs, t2);
    if s2.abs() > 1e-12 {
        s1 / s2
    } else {
        1.0
    }
}

/// Estimate LAXPC30 gas density fraction `rho` (relative to launch) at `ndays` (`backshiftv3.f:428-430`).
pub fn estimate_lxp30_density(ndays: f64) -> f64 {
    let mut rho = 1.0;
    if ndays > 130.0 {
        rho = 1.0 - (ndays - 130.0) * 1.50e-3;
    }
    if ndays > 700.0 {
        rho = 0.4 * (-(ndays - 500.0) / 200.0).exp();
    }
    rho.max(0.01)
}

/// Recommend the best matching TIFR CALDB `.rmf` filename for a given unit, observation day,
/// 30 keV peak channel `chs`, layer selection, and single-event mode (`backshiftv3.f:417-526`).
pub fn recommend_rmf(
    unit: LaxpcUnit,
    ndays: f64,
    source_peak_ch: f64,
    layer1_only: bool,
    single_only: bool,
) -> String {
    let gar = unit.reference_30kev_channel(ndays);
    let mut dch = source_peak_ch - gar;

    // LAXPC30 density rounding to steps of 5% (`backshiftv3.f:434-440`)
    let rho = estimate_lxp30_density(ndays);
    let mut irho = (((rho * 100.0 + 2.5) as i32) / 5) * 5;
    if irho > 99 {
        irho = 100;
    }
    if rho < 0.08 {
        irho = 6;
    }
    if rho < 0.05 {
        irho = 4;
    }

    let mut idel = (dch.abs() + 0.5) as i32;
    let mut id2 = (((dch.abs() + 1.0) / 2.0) as i32) * 2;
    let mut id4 = (((dch.abs() + 2.0) / 4.0) as i32) * 4;
    let mut id5 = (((dch.abs() + 2.5) / 5.0) as i32) * 5;

    if unit != LaxpcUnit::Lxp30 && idel > 20 {
        idel = 20;
    }
    if unit == LaxpcUnit::Lxp30 {
        if id2 > 30 {
            id2 = 30;
        }
        if id4 > 32 {
            id4 = 32;
        }
        if id5 > 95 {
            id5 = 95;
        }
    }

    // Post-anomaly LAXPC10 large negative shift offsets (`backshiftv3.f:452-467`)
    if unit == LaxpcUnit::Lxp10 && ndays > 917.0 {
        let offset = if ndays <= 991.0 {
            140.0
        } else if ndays <= 1036.0 {
            175.0
        } else if ndays <= 1112.0 {
            197.0
        } else if ndays <= 2944.0 {
            212.0
        } else {
            220.0
        };
        dch -= offset;
        idel = (dch.abs() + 0.5) as i32;
    }

    let prefix = unit.short_name();
    let layer_tag = if layer1_only { "L1" } else { "" };
    let se_tag = if single_only { "SE" } else { "" };
    let rho_tag = if unit == LaxpcUnit::Lxp30 && irho < 97 {
        format!("rho{irho:02}")
    } else {
        String::new()
    };

    if idel == 0 || (unit == LaxpcUnit::Lxp30 && dch.abs() < 1.0) {
        let ver = if unit == LaxpcUnit::Lxp10 { "v1.1" } else { "v1.0" };
        return format!("{prefix}{rho_tag}{layer_tag}{se_tag}{ver}.rmf");
    }

    if dch < 0.0 {
        if idel >= 100 {
            return format!("{prefix}cshm{idel:03}{layer_tag}{se_tag}v2.2.rmf");
        }
        let step = if unit == LaxpcUnit::Lxp30 {
            if irho < 84 { id4 } else { id2 }
        } else {
            idel
        };
        let ver = if unit == LaxpcUnit::Lxp10 { "v1.1" } else { "v1.0" };
        format!("{prefix}cshm{step:02}{rho_tag}{layer_tag}{se_tag}{ver}.rmf")
    } else {
        let step = if unit == LaxpcUnit::Lxp30 {
            if irho < 84 {
                id4
            } else if dch > 32.5 {
                id5
            } else {
                id2
            }
        } else {
            idel
        };
        let ver = if unit == LaxpcUnit::Lxp10 { "v1.1" } else { "v1.0" };
        format!("{prefix}cshp{step:02}{rho_tag}{layer_tag}{se_tag}{ver}.rmf")
    }
}

/// Compute the full gain-shift solution for aligning a background spectrum from `background_day`
/// (with known 30 keV peak `background_peak_ch`, or interpolated from `gain_table` if `None`)
/// to a source observation at `source_day`.
#[allow(clippy::too_many_arguments)]
pub fn compute_gain_shift(
    unit: LaxpcUnit,
    source_day: f64,
    background_day: f64,
    background_peak_ch: Option<f64>,
    gain_table: &GainTable,
    hv_table: Option<&HvJumpTable>,
    layer1_only: bool,
    single_only: bool,
) -> Result<GainShiftSolution> {
    let mut source_peak_ch = gain_table.peak_30kev_at(unit, source_day, hv_table)?;
    let mut chb = match background_peak_ch {
        Some(c) => c,
        None => gain_table.peak_30kev_at(unit, background_day, hv_table)?,
    };

    // Match `backshiftv3.f:298-307`: apply LAXPC20 linear drift fit ONLY when BOTH
    // `source_day` and `background_day` exceed 550 (or both exceed 2977).
    if unit == LaxpcUnit::Lxp20 && source_day > 2977.0 && background_day > 2977.0 {
        source_peak_ch = 228.2 - 0.00571 * (source_day - 2977.0);
        chb = 228.2 - 0.00571 * (background_day - 2977.0);
    } else if unit == LaxpcUnit::Lxp20 && source_day > 550.0 && background_day > 550.0 {
        source_peak_ch = 230.1 - 0.00571 * source_day;
        chb = 230.1 - 0.00571 * background_day;
    }

    let delta_ch_30kev = source_peak_ch - chb;
    let gain_ratio = if chb.abs() > 1e-12 {
        source_peak_ch / chb
    } else {
        1.0
    };
    let secular_scale = secular_background_scale(unit, source_day, background_day);
    let crossed_hv_jumps = hv_table
        .map(|hv| hv.jumps_between(unit, background_day, source_day))
        .unwrap_or_default();
    let lxp20_gain_offset_kev = if unit == LaxpcUnit::Lxp20 {
        Some(0.08919 - 0.0003001 * source_day)
    } else {
        None
    };
    let recommended_rmf = recommend_rmf(unit, source_day, source_peak_ch, layer1_only, single_only);

    Ok(GainShiftSolution {
        unit,
        source_day,
        background_day,
        source_peak_ch,
        background_peak_ch: chb,
        delta_ch_30kev,
        gain_ratio,
        secular_scale,
        crossed_hv_jumps,
        lxp20_gain_offset_kev,
        recommended_rmf,
    })
}

/// Shift a multi-anode background spectrum (`[nch][11]` rates: A1..A10 + sum A1..A7) by `delta_ch_30kev`
/// and apply `secular_scale`, matching `backshiftv3.f:gainsh` (lines 661-726).
pub fn shift_background_spectrum(
    unit: LaxpcUnit,
    delta_ch_30kev: f64,
    secular_scale: f64,
    ndays: f64,
    channel_map: &ChannelMap,
    spec: &[[f64; 11]],
) -> Result<Vec<[f64; 11]>> {
    let nch = spec.len();
    if nch == 0 || channel_map.len() < nch {
        return Err(LaxpcError::Background(format!(
            "spectrum length ({nch}) incompatible with channel map length ({})",
            channel_map.len()
        )));
    }

    let nf = (1024 / nch).max(1) as f64;
    let e2 = unit.gain_shift_e2();
    let mut icu = (unit.uld_cutoff_channel_1024() + 1) / (nf as usize);
    if unit == LaxpcUnit::Lxp10 && ndays > 917.0 {
        icu = 158;
    }
    icu = icu.min(nch);

    let de = (delta_ch_30kev / 30.0) / nf;
    let mut shifted = vec![[0.0_f64; 11]; nch];
    let mut last_n1 = 1usize;

    // 1-based channel loop `i = 1..=icu` matching `backshiftv3.f:693-712`
    for i1 in 1..=icu {
        let el = channel_map.e_min_kev[i1 - 1];
        let eu = channel_map.e_max_kev[i1 - 1];
        let mut dn1 = de * el * (1.0 + e2 * el);
        let dn2 = 1.0 + de * eu * (1.0 + e2 * eu);
        let dn = dn2 - dn1;
        if dn.abs() <= 1e-15 {
            continue;
        }
        let n1 = ((i1 as f64 + dn1).floor() as isize).max(1) as usize;
        let n2 = (((i1 as f64 + 1.0 + dn2).floor() as isize).max(icu as isize) as usize).min(nch);
        last_n1 = n1;

        for j1 in n1..=n2 {
            let dj = ((j1 as f64 + 1.0) - i1 as f64).min(dn2);
            if dj > dn1 {
                let weight = (dj - dn1) / dn;
                for k in 0..11 {
                    shifted[j1 - 1][k] += weight * spec[i1 - 1][k];
                }
            }
            dn1 = (j1 as f64 + 1.0) - i1 as f64;
        }
    }

    // Fill edge gap below ULD when shift is negative (`backshiftv3.f:714-720`)
    if de < 0.0 && last_n1 >= 2 {
        for i1 in (last_n1 - 1).max(2)..=icu {
            for k in 0..11 {
                shifted[i1 - 1][k] = shifted[i1 - 2][k];
            }
        }
    }

    // Copy channels above ULD unmodified (`backshiftv3.f:721-725`)
    if icu < nch {
        shifted[icu..nch].copy_from_slice(&spec[icu..nch]);
    }

    // Apply secular scaling and recompute main-anode sum column 10 (A1-A7)
    for row in &mut shifted {
        let mut sum_main = 0.0;
        for (k, cell) in row.iter_mut().take(10).enumerate() {
            *cell *= secular_scale;
            if k < 7 {
                sum_main += *cell;
            }
        }
        row[10] = sum_main;
    }

    Ok(shifted)
}

/// Compute the fractional background light-curve gain-shift adjustment `fshbk = (c_shifted - c_orig) / c_orig`
/// for a chosen anode selection and native channel range `[ch_min, ch_max]`.
pub fn background_lightcurve_shift_fraction(
    orig: &[[f64; 11]],
    shifted: &[[f64; 11]],
    anodes: &AnodeSelection,
    ch_min: usize,
    ch_max: usize,
) -> (f64, f64) {
    let mut c_orig = 0.0;
    let mut c_shifted = 0.0;
    let n = orig.len().min(shifted.len());
    let hi = ch_max.min(n.saturating_sub(1));

    for i in ch_min..=hi {
        for anode in 1..=10u8 {
            if anodes.matches(anode) {
                let k = (anode - 1) as usize;
                c_orig += orig[i][k];
                c_shifted += shifted[i][k];
            }
        }
    }

    let delta = c_shifted - c_orig;
    let frac = if c_orig > 0.0 { delta / c_orig } else { 0.0 };
    (delta, frac)
}

/// Shift a background light curve in place by applying the gain-shift factor `1.0 + fshbk`
/// (`backshiftv3.f:606`).
pub fn apply_lightcurve_gain_shift(bkg_rates: &mut [f64], fshbk: f64) {
    let mult = 1.0 + fshbk;
    for r in bkg_rates.iter_mut() {
        *r *= mult;
    }
}

/// Parsed ASCII background spectrum (`backlxp{10,20,30}*.pha`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackgroundPha {
    pub exposure: f64,
    pub mission_day: f64,
    /// Per-channel rates `[nch][11]` (columns 0..10 = A1..A10, column 10 = sum A1..A7).
    pub rates: Vec<[f64; 11]>,
}

impl BackgroundPha {
    /// Read a TIFR `.pha` ASCII background spectrum file (`backlxp?0*.pha`).
    pub fn from_file(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let header = lines
            .next()
            .ok_or_else(|| LaxpcError::Background("empty background .pha file".to_string()))??;

        // Example: `#  Exposure =   1.204432E+04  lxp1  Channels:    0  511  day:  100  rate: ...`
        let exposure = parse_token_after(&header, "Exposure =")
            .or_else(|| parse_token_after(&header, "Exposure="))
            .unwrap_or(1.0);
        let mission_day = parse_token_after(&header, "day:").unwrap_or(0.0);

        let mut rates = Vec::with_capacity(512);
        for line_res in lines {
            let line = line_res?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let vals: Vec<f64> = trimmed
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if vals.len() >= 12 {
                let mut row = [0.0_f64; 11];
                row.copy_from_slice(&vals[1..12]);
                rates.push(row);
            }
        }
        if rates.is_empty() {
            return Err(LaxpcError::Background(format!(
                "no spectral rows in {}",
                path.display()
            )));
        }
        Ok(Self {
            exposure,
            mission_day,
            rates,
        })
    }
}

fn parse_token_after(haystack: &str, needle: &str) -> Option<f64> {
    let idx = haystack.find(needle)?;
    let rest = &haystack[idx + needle.len()..];
    rest.split_whitespace().next()?.parse::<f64>().ok()
}

/// 2D Earth longitude/latitude B-spline + 5-harmonic diurnal background model
/// parsed from Fortran unformatted `backfit.lxp?0*` files (`laxpcl1.f:1160-1170,1433-1447`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackgroundModel2D {
    /// Longitude knots (`no`, typically 40, spanning 0..360 deg).
    pub lon_knots: Vec<f64>,
    /// Latitude knots (`na`, typically 10, spanning -6..+6 deg).
    pub lat_knots: Vec<f64>,
    /// Leading dimension `iw` (`no + 2`) of the 2D B-spline coefficient array.
    pub iw: usize,
    /// Full coefficient vector: first `(no+2)*(na+2)` are column-major 2D B-spline weights,
    /// followed by `nv` (typically 10 = 5 harmonics x (sin, cos)) diurnal harmonic coefficients.
    pub coeffs: Vec<f64>,
    /// Mean background count rate `cav` (cps) used for +/- 40 cps clamping.
    pub mean_rate: f64,
    /// Diurnal oscillation period `per` (s, typically ~84494.9 s).
    pub diurnal_period_sec: f64,
}

impl BackgroundModel2D {
    /// Load a binary Fortran unformatted `backfit.*` file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let mut file = File::open(path)?;
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes)?;
        Self::from_fortran_bytes(&bytes)
    }

    /// Parse Fortran unformatted sequential records (4-byte LE length marker before & after each record).
    pub fn from_fortran_bytes(bytes: &[u8]) -> Result<Self> {
        let records = read_fortran_records(bytes)?;
        if records.len() < 4 {
            return Err(LaxpcError::Background(format!(
                "backfit file expected >= 4 Fortran records, found {}",
                records.len()
            )));
        }

        // Record 0: `no, na, iw, mbk` (4 x i32) + `cav, per` (2 x f64) [+ optional trailing metadata]
        let r0 = records[0];
        if r0.len() < 32 {
            return Err(LaxpcError::Background(format!(
                "backfit record 0 too short: {} bytes (expected >= 32)",
                r0.len()
            )));
        }
        let no_raw = i32::from_le_bytes(r0[0..4].try_into().unwrap());
        let na_raw = i32::from_le_bytes(r0[4..8].try_into().unwrap());
        let iw_raw = i32::from_le_bytes(r0[8..12].try_into().unwrap());
        let mbk_raw = i32::from_le_bytes(r0[12..16].try_into().unwrap());
        let mean_rate = f64::from_le_bytes(r0[16..24].try_into().unwrap());
        let diurnal_period_sec = f64::from_le_bytes(r0[24..32].try_into().unwrap());

        if no_raw < 2 || na_raw < 2 || iw_raw <= 0 || mbk_raw <= 0 {
            return Err(LaxpcError::Background(format!(
                "invalid backfit dimensions (negative or too small): no={no_raw}, na={na_raw}, iw={iw_raw}, mbk={mbk_raw}"
            )));
        }
        let no = no_raw as usize;
        let na = na_raw as usize;
        let iw = iw_raw as usize;
        let mbk = mbk_raw as usize;

        let min_spline = no
            .checked_add(2)
            .and_then(|a| na.checked_add(2).and_then(|b| a.checked_mul(b)))
            .ok_or_else(|| LaxpcError::Background("backfit dimension overflow".to_string()))?;

        if no > 10_000 || na > 10_000 || mbk > 10_000_000 || mbk < min_spline {
            return Err(LaxpcError::Background(format!(
                "invalid backfit dimensions: no={no}, na={na}, iw={iw}, mbk={mbk}"
            )));
        }

        let lon_knots = read_f64_vec(records[1], no)?;
        let lat_knots = read_f64_vec(records[2], na)?;
        let coeffs = read_f64_vec(records[3], mbk)?;

        Ok(Self {
            lon_knots,
            lat_knots,
            iw,
            coeffs,
            mean_rate,
            diurnal_period_sec: if diurnal_period_sec > 0.0 {
                diurnal_period_sec
            } else {
                DEFAULT_DIURNAL_PERIOD_SEC
            },
        })
    }

    /// Evaluate the 2D B-spline + 5-harmonic diurnal background count rate (cps) at spacecraft
    /// `(lon_deg, lat_deg)` and AstroSat UTC time `t_utc` (`laxpcl1.f:1433-1447`).
    pub fn evaluate(&self, lon_deg: f64, lat_deg: f64, t_utc: f64) -> f64 {
        let mut lon = lon_deg;
        if lon < 0.0 {
            lon += 360.0;
        }
        if lon > 252.0 {
            return (self.mean_rate - 40.0).max(0.0);
        }

        let no = self.lon_knots.len();
        let na = self.lat_knots.len();
        let mut bks = eval_bspline_2d(
            &self.lon_knots,
            &self.lat_knots,
            4,
            &self.coeffs,
            self.iw,
            lon,
            lat_deg,
        );

        // Add diurnal harmonics (`laxpcl1.f:1439-1443`):
        // `nv = mbk - (no+2)*(na+2)`, pairs `(sin(h*om*tu), cos(h*om*tu))` for h = 1..=nv/2
        let spline_count = (no + 2) * (na + 2);
        if self.coeffs.len() > spline_count {
            let nv = self.coeffs.len() - spline_count;
            let om = 2.0 * std::f64::consts::PI / self.diurnal_period_sec;
            let tu = t_utc.rem_euclid(self.diurnal_period_sec);
            let mut j = 0usize;
            while j + 1 < nv {
                let harmonic = (j / 2 + 1) as f64;
                let phase = harmonic * om * tu;
                bks += self.coeffs[spline_count + j] * phase.sin()
                    + self.coeffs[spline_count + j + 1] * phase.cos();
                j += 2;
            }
        }

        bks.clamp(self.mean_rate - 40.0, self.mean_rate + 40.0)
    }
}

/// Evaluate cubic (or arbitrary order `k`) B-spline basis functions at `x0` over knot vector `knots`,
/// matching `laxpcl1.f:BSPLIN` (lines 3436-3581).
///
/// Returns `(left_zero_based, basis_k)` where `basis_k[r]` is $B_{\text{left}+r, k}(x_0)$ for `r in 0..k`.
pub fn eval_bspline_basis(knots: &[f64], k: usize, x0: f64) -> (usize, Vec<f64>) {
    let nx = knots.len();
    if nx < 2 || k == 0 {
        return (0, vec![0.0; k]);
    }

    // Bracket `x0` in `knots` (0-based `low` in `0..=nx-2`)
    let low = if x0 <= knots[0] {
        0
    } else if x0 >= knots[nx - 1] {
        nx - 2
    } else {
        let mut lo = 0usize;
        let mut hi = nx - 1;
        while hi - lo > 1 && x0 != knots[lo] {
            let mid = (lo + hi) / 2;
            if x0 >= knots[mid] {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    };

    // Extended clamped knot vector `t`: indices `-k..nx+k` shifted by `+k`
    let knot_at = |idx: isize| -> f64 {
        if idx < 0 {
            knots[0]
        } else if (idx as usize) >= nx {
            knots[nx - 1]
        } else {
            knots[idx as usize]
        }
    };

    let mut b = vec![0.0_f64; k];
    let mut dr = vec![0.0_f64; k];
    let mut dl = vec![0.0_f64; k];
    b[0] = 1.0;

    let low_i = low as isize;
    for j in 1..k {
        dr[j - 1] = knot_at(low_i + j as isize) - x0;
        dl[j - 1] = x0 - knot_at(low_i + 1 - j as isize);
        let mut saved = 0.0;
        for r in 0..j {
            let denom = dr[r] + dl[j - 1 - r];
            let temp = if denom.abs() > 0.0 { b[r] / denom } else { 0.0 };
            b[r] = saved + temp * dr[r];
            saved = temp * dl[j - 1 - r];
        }
        b[j] = saved;
    }

    (low, b)
}

/// Evaluate a 2D tensor-product B-spline expansion `sum WT(i,j) phi_i(x0) psi_j(y0)`
/// stored in Fortran column-major order with leading dimension `iw`, matching `laxpcl1.f:BSPEV2`.
pub fn eval_bspline_2d(
    x_knots: &[f64],
    y_knots: &[f64],
    k: usize,
    wt: &[f64],
    iw: usize,
    x0: f64,
    y0: f64,
) -> f64 {
    let (lx, bx) = eval_bspline_basis(x_knots, k, x0);
    let (ly, by) = eval_bspline_basis(y_knots, k, y0);

    let mut sum = 0.0;
    for (rx, &bx_val) in bx.iter().take(k).enumerate() {
        let ix = lx + rx;
        for (ry, &by_val) in by.iter().take(k).enumerate() {
            let iy = ly + ry;
            let flat = iy * iw + ix;
            if flat < wt.len() {
                sum += wt[flat] * bx_val * by_val;
            }
        }
    }
    sum
}

/// Parse Fortran sequential unformatted records `[len: u32 LE][payload: len bytes][len: u32 LE]`.
fn read_fortran_records(bytes: &[u8]) -> Result<Vec<&[u8]>> {
    let mut out = Vec::new();
    let mut pos = 0usize;
    while pos + 4 <= bytes.len() {
        let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        if pos + len + 4 > bytes.len() {
            return Err(LaxpcError::Background(format!(
                "truncated Fortran record at offset {pos}: declared {len} bytes"
            )));
        }
        let payload = &bytes[pos..pos + len];
        pos += len;
        let tail = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        if tail != len {
            return Err(LaxpcError::Background(format!(
                "mismatched Fortran record trailer: header={len}, trailer={tail}"
            )));
        }
        out.push(payload);
    }
    Ok(out)
}

fn read_f64_vec(rec: &[u8], count: usize) -> Result<Vec<f64>> {
    let needed = count
        .checked_mul(8)
        .ok_or_else(|| LaxpcError::Background(format!("f64 count overflow: {count}")))?;
    if rec.len() < needed {
        return Err(LaxpcError::Background(format!(
            "record too short for {count} f64 values: {} bytes (need {needed})",
            rec.len()
        )));
    }
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        out.push(f64::from_le_bytes(rec[i * 8..(i + 1) * 8].try_into().unwrap()));
    }
    Ok(out)
}
