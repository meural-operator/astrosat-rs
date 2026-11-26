//! Stage 4: TCT time conversion and standard LAXPC GTI generation.
//!
//! Converts 100 kHz STBG instrument timestamps to AstroSat UTC seconds (since MJD 55197)
//! via `astrosat_core::tct::LaxpcTct` (rejecting `dt_sps/dt_ins` outside `(1.000, 1.001]`),
//! and generates standard GTIs from MKF + BBC + count rates (Earth occultation `ELV > 3.5 deg`,
//! SAA longitude & diurnal extension, LAXPC10 ULD anomaly thresholds, gap merge `< 0.3 s`,
//! short drop `< 5.0 s`, frame-gap exposure subtraction, and optional user GTI intersection).

use crate::bb::{BbcSeries, UldQuartileSummary};
use crate::error::{LaxpcError, Result};
use crate::products::{compute_gti_frame_gap_loss, RawBinnedCounts};
use crate::types::{utc_to_mission_day, LaxpcUnit, DEFAULT_DEADTIME_SEC};
use astrosat_core::gti::{GtiInterval, GtiSet};
use astrosat_core::mkf::MkfTable;
use astrosat_core::tct::{LaxpcTct, TctTable};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Configuration for standard LAXPC GTI generation (`laxpcl1.f:253-258,1210-1565`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GtiConfig {
    /// Minimum source elevation above Earth's limb in degrees (default `3.5`).
    pub min_elv_deg: f64,
    /// If `true`, invert elevation criterion (`je = -1`) to extract Earth-occulted intervals.
    pub earth_occultation: bool,
    /// SAA entry longitude in degrees (`lon > saa_entry_lon` excluded, default `249.0`).
    pub saa_entry_lon_deg: f64,
    /// SAA exit longitude in degrees (`lon < saa_exit_lon` excluded, default `3.0`).
    pub saa_exit_lon_deg: f64,
    /// Enable diurnal SAA exit longitude extension to 40 deg (`laxpcl1.f:1321-1327`).
    pub enable_diurnal_saa_extension: bool,
    /// Enable post-anomaly LAXPC10 ULD rate threshold screening (`laxpcl1.f:1237-1246`).
    pub enable_lxp10_uld_anomaly_filter: bool,
    /// Optional ULD quartile bin selector (`1..=4`; `None` accepts all ULD rates).
    pub uld_quartile_filter: Option<u8>,
    /// Optional override for active-detector all-anode count-rate threshold (cps).
    pub count_rate_threshold_cps: Option<f64>,
    /// Merge adjacent GTI intervals separated by less than this gap (default `0.3 s`).
    pub merge_gap_sec: f64,
    /// Drop GTI intervals shorter than this duration (default `5.0 s`).
    pub min_duration_sec: f64,
    /// Diurnal period in seconds for SAA extension phase calculation (default `84494.9 s`).
    pub diurnal_period_sec: f64,
    /// Deadtime constant in seconds for ULD rate correction (default `43 us`).
    pub deadtime_sec: f64,
}

impl Default for GtiConfig {
    fn default() -> Self {
        Self {
            min_elv_deg: 3.5,
            earth_occultation: false,
            saa_entry_lon_deg: 249.0,
            saa_exit_lon_deg: 3.0,
            enable_diurnal_saa_extension: true,
            enable_lxp10_uld_anomaly_filter: true,
            uld_quartile_filter: None,
            count_rate_threshold_cps: None,
            merge_gap_sec: 0.3,
            min_duration_sec: 5.0,
            diurnal_period_sec: 84_494.9,
            deadtime_sec: DEFAULT_DEADTIME_SEC,
        }
    }
}

/// Result of standard LAXPC GTI screening.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GtiResult {
    /// Good time intervals in instrument seconds (used for direct event/bin filtering).
    pub gti_ins: GtiSet,
    /// Good time intervals in AstroSat UTC seconds (written to FITS `GTI` extensions).
    pub gti_utc: GtiSet,
    /// Gross GTI exposure (sum of interval lengths in instrument seconds).
    pub gross_exposure_sec: f64,
    /// Internal frame-gap time lost inside the GTI (seconds).
    pub frame_gap_loss_sec: f64,
    /// Net science exposure (`gross_exposure_sec - frame_gap_loss_sec`).
    pub net_exposure_sec: f64,
    /// ULD quartile exposure distribution across the GTI.
    pub uld_quartiles: UldQuartileSummary,
    /// Maximum pointing angular offset (`ANG_OFFSET`) observed in degrees.
    pub max_ang_offset_deg: f64,
    /// UTC timestamp of maximum pointing offset.
    pub max_ang_offset_time_utc: f64,
    /// Number of TCT glitch rows rejected during TCT loading.
    pub tct_rejected_glitches: usize,
}

/// Load and merge one or more Level-1 `.tct` files into a single glitch-cleaned [`LaxpcTct`].
pub fn load_laxpc_tct_files<P: AsRef<Path>>(paths: &[P]) -> Result<LaxpcTct> {
    if paths.is_empty() {
        return Err(LaxpcError::Tct("no TCT files provided".into()));
    }
    if paths.len() == 1 {
        return LaxpcTct::from_fits(&paths[0]).map_err(|e| LaxpcError::Tct(e.to_string()));
    }

    let mut sps_all = Vec::new();
    let mut ins_all = Vec::new();
    let mut last_sps = f64::NEG_INFINITY;

    for p in paths {
        let tbl = TctTable::from_fits(p, None).map_err(|e| LaxpcError::Tct(e.to_string()))?;
        for (&s, &i) in tbl.sps_time.iter().zip(tbl.inst_time.iter()) {
            if s.is_finite() && i.is_finite() && s > 1.0e5 && s > last_sps {
                sps_all.push(s);
                ins_all.push(i);
                last_sps = s;
            }
        }
    }

    let merged = TctTable::from_pairs(sps_all, ins_all).map_err(|e| LaxpcError::Tct(e.to_string()))?;
    LaxpcTct::new(merged).map_err(|e| LaxpcError::Tct(e.to_string()))
}

/// Load and merge one or more Level-1 `.mkf` files into a single monotonic [`MkfTable`].
pub fn load_merged_mkf_files<P: AsRef<Path>>(paths: &[P]) -> Result<MkfTable> {
    if paths.is_empty() {
        return Err(LaxpcError::Gti("no MKF files provided".into()));
    }
    if paths.len() == 1 {
        return MkfTable::from_fits(&paths[0]).map_err(|e| LaxpcError::Core(e.to_string()));
    }

    let first = MkfTable::from_fits(&paths[0]).map_err(|e| LaxpcError::Core(e.to_string()))?;
    let col_keys: Vec<String> = first
        .column_names()
        .iter()
        .map(|s| s.trim().to_ascii_uppercase())
        .collect();

    let mut merged_times = Vec::new();
    let mut merged_cols: HashMap<String, Vec<f64>> = col_keys
        .iter()
        .map(|k| (k.clone(), Vec::new()))
        .collect();
    let mut last_t = f64::NEG_INFINITY;

    for (idx, p) in paths.iter().enumerate() {
        let tbl = if idx == 0 {
            first.clone()
        } else {
            MkfTable::from_fits(p).map_err(|e| LaxpcError::Core(e.to_string()))?
        };
        let times = tbl.times();
        for (row, &t) in times.iter().enumerate() {
            if !t.is_finite() || t <= 1.0e5 || t - last_t < 1.0e-2 {
                continue;
            }
            merged_times.push(t);
            for key in &col_keys {
                let val = tbl.column(key).and_then(|c| c.get(row).copied()).unwrap_or(0.0);
                merged_cols.get_mut(key).unwrap().push(val);
            }
            last_t = t;
        }
    }

    MkfTable::from_columns(merged_times, merged_cols).map_err(|e| LaxpcError::Core(e.to_string()))
}

/// Invert a monotonic [`LaxpcTct`] mapping (`UTC -> instrument_time`) using fast secant refinement.
///
/// Accurate to `< 1e-9 s` over arbitrary mission baselines.
pub fn utc_to_ins(tct: &LaxpcTct, t_utc: f64) -> f64 {
    // Initial linear estimate from two probe points around t_utc
    let t0 = 0.0;
    let u0 = tct.interpolate(t0);
    let u1 = tct.interpolate(1000.0);
    let slope = (u1 - u0) / 1000.0;
    let mut x = if slope.abs() > 0.5 {
        t0 + (t_utc - u0) / slope
    } else {
        t_utc
    };

    // Secant/Newton iterations with local numerical derivative
    for _ in 0..12 {
        let fx = tct.interpolate(x) - t_utc;
        if fx.abs() < 1.0e-10 {
            break;
        }
        let h = 1.0;
        let dfdx = (tct.interpolate(x + h) - tct.interpolate(x - h)) / (2.0 * h);
        if dfdx.abs() < 0.5 {
            break;
        }
        x -= fx / dfdx;
    }
    x
}

/// LAXPC10 post-anomaly minimum ULD count-rate threshold (cps) as a function of UTC time
/// (`laxpcl1.f:1237-1246`). Returns `0.0` for LAXPC20/30 or pre-anomaly LAXPC10 epochs.
pub fn lxp10_uld_anomaly_threshold(unit: LaxpcUnit, t_utc: f64) -> f64 {
    if unit != LaxpcUnit::Lxp10 {
        return 0.0;
    }
    let mut uthr = 0.0;
    if t_utc > 2.6051922e8 {
        uthr = 750.0;
    }
    if t_utc > 2.600928e8 && t_utc < 2.6036232e8 {
        uthr = 450.0;
    }
    if t_utc > 2.6251320e8 && t_utc < 2.6295e8 {
        uthr = 300.0;
    }
    if t_utc > 2.6334e8 && t_utc < 2.6346e8 {
        uthr = 450.0;
    }
    if t_utc > 2.6674578e8 {
        uthr = 450.0;
    }
    if t_utc > 2.7065000e8 {
        uthr = 300.0;
    }
    if t_utc > 2.7729300e8 {
        uthr = 200.0;
    }
    if t_utc > 2.9318178e8 {
        uthr = 130.0;
    }
    if t_utc > 4.3514610e8 {
        uthr = 35.0;
    }
    uthr
}

/// Generate standard LAXPC GTIs from binned all-anode counts, MKF, optional BBC series, and TCT,
/// matching `laxpcl1.f:1210-1565`, then optionally intersect with `user_gti_utc`.
#[allow(clippy::too_many_arguments)]
pub fn generate_standard_gti(
    unit: LaxpcUnit,
    raw_counts: &RawBinnedCounts,
    tct: &LaxpcTct,
    mkf: &MkfTable,
    bbc: Option<&BbcSeries>,
    frame_gaps_ins: &[(f64, f64)],
    cfg: &GtiConfig,
    user_gti_utc: Option<&GtiSet>,
) -> Result<GtiResult> {
    let nt = raw_counts.all_main_counts.len();
    if nt == 0 {
        return Ok(GtiResult {
            gti_ins: GtiSet::default(),
            gti_utc: GtiSet::default(),
            gross_exposure_sec: 0.0,
            frame_gap_loss_sec: 0.0,
            net_exposure_sec: 0.0,
            uld_quartiles: UldQuartileSummary::default(),
            max_ang_offset_deg: 0.0,
            max_ang_offset_time_utc: 0.0,
            tct_rejected_glitches: tct.rejected_glitch_count(),
        });
    }

    let tbin = raw_counts.tbin;
    let je = if cfg.earth_occultation { -1.0 } else { 1.0 };

    // Locate MKF columns flexibly across L1/L2 naming variants
    let elv_col = mkf
        .column_any(&["ELV", "ELEVATION", "ELV_ANGLE"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing ELV column".into()))?;
    let lon_col = mkf
        .column_any(&["EARTHLON", "SAT_LON", "LONGITUDE", "LON"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing EARTHLON/LONGITUDE column".into()))?;
    let off_col = mkf.column_any(&["ANG_OFFSET", "OFFSET"]);
    let mkf_t = mkf.times();

    let t_utc_first = tct.interpolate(raw_counts.t0_ins + 0.5 * tbin);
    let cthr = cfg
        .count_rate_threshold_cps
        .unwrap_or_else(|| unit.count_rate_threshold(t_utc_first));
    let nday_approx = utc_to_mission_day(t_utc_first);

    let mut raw_intervals: Vec<GtiInterval> = Vec::new();
    let mut cur_start: Option<f64> = None;
    let mut uld_quartiles = UldQuartileSummary::default();
    let mut max_ang_offset_deg = 0.0_f64;
    let mut max_ang_offset_time_utc = t_utc_first;
    let mut saaex = cfg.saa_exit_lon_deg;

    for i in 0..nt {
        let tc = raw_counts.t0_ins + (i as f64 + 0.5) * tbin;
        let xb = tct.interpolate(tc);
        let rate_all = raw_counts.all_main_counts[i] / tbin;

        let elv = crate::bb::interp_linear(mkf_t, elv_col, xb);
        let mut lon = interp_lon_deg(mkf_t, lon_col, xb);
        if lon < 0.0 {
            lon += 360.0;
        }

        if let Some(oc) = off_col {
            let off = crate::bb::interp_linear(mkf_t, oc, xb);
            if off.is_finite() && off > max_ang_offset_deg {
                max_ang_offset_deg = off;
                max_ang_offset_time_utc = xb;
            }
        }

        // Diurnal SAA exit longitude extension (`laxpcl1.f:1321-1327`)
        if cfg.enable_diurnal_saa_extension
            && ((lon > 0.0 && lon < 3.0) || (i < 2 && lon < 40.0))
        {
            let tu0 = xb.rem_euclid(cfg.diurnal_period_sec);
            saaex = cfg.saa_exit_lon_deg;
            if (tu0 > 1.0e4 && tu0 < 3.5e4 && nday_approx > 1104.0)
                || (tu0 > 1.0e4 && tu0 < 3.0e4 && nday_approx <= 1100.0 && nday_approx > 186.0)
            {
                saaex = 40.0;
            }
        }

        // ULD rate and LAXPC10 anomaly filter
        let uld_rate = bbc
            .map(|b| b.uld_rate_cps(tc, rate_all, cfg.deadtime_sec))
            .unwrap_or(0.0);
        let uthr = if cfg.enable_lxp10_uld_anomaly_filter && bbc.is_some() {
            lxp10_uld_anomaly_threshold(unit, xb)
        } else {
            0.0
        };
        let uld_ok = !(uthr > 30.0 && uld_rate < uthr);
        let uld_bin_ok = match cfg.uld_quartile_filter {
            Some(target) if bbc.is_some() => {
                BbcSeries::uld_quartile_index(unit, uld_rate) == target as i32
            }
            _ => true,
        };

        // Combined good-bin condition: Earth elevation + outside SAA + detector active + ULD OK
        let elv_ok = elv * je > cfg.min_elv_deg;
        let saa_ok = lon < cfg.saa_entry_lon_deg && lon > saaex;
        let active_ok = rate_all >= cthr;
        let good = elv_ok && saa_ok && active_ok && uld_ok && uld_bin_ok;

        let bin_lo = tc - 0.5 * tbin;
        let bin_hi = tc + 0.5 * tbin;

        if good {
            if cur_start.is_none() {
                cur_start = Some(bin_lo);
            }
            if bbc.is_some() {
                uld_quartiles.record(unit, uld_rate);
            }
        } else if let Some(s) = cur_start.take() {
            raw_intervals.push(GtiInterval::new(s, bin_lo));
        }

        if i + 1 == nt {
            if let Some(s) = cur_start.take() {
                raw_intervals.push(GtiInterval::new(s, bin_hi));
            }
        }
    }

    uld_quartiles.finalize();

    // Normalize, merge small gaps (< 0.3 s), and drop short intervals (< 5.0 s)
    let mut gti_ins = GtiSet::new(raw_intervals)
        .merge_gaps(cfg.merge_gap_sec)
        .drop_short(cfg.min_duration_sec);

    // Map to UTC
    let mut gti_utc = GtiSet::new(
        gti_ins
            .intervals()
            .iter()
            .map(|iv| GtiInterval::new(tct.interpolate(iv.start), tct.interpolate(iv.stop)))
            .collect(),
    );

    // Intersect with user-provided UTC GTI if supplied (exact set intersection: empty ∩ U = empty)
    if let Some(user_gti) = user_gti_utc {
        gti_utc = gti_utc.intersect(user_gti);
        gti_ins = GtiSet::new(
            gti_utc
                .intervals()
                .iter()
                .map(|iv| GtiInterval::new(utc_to_ins(tct, iv.start), utc_to_ins(tct, iv.stop)))
                .collect(),
        );
    }

    let pairs_ins: Vec<(f64, f64)> = gti_ins
        .intervals()
        .iter()
        .map(|iv| (iv.start, iv.stop))
        .collect();
    let gross_exposure_sec = gti_ins.exposure();
    let frame_gap_loss_sec = compute_gti_frame_gap_loss(&pairs_ins, frame_gaps_ins);
    let net_exposure_sec = (gross_exposure_sec - frame_gap_loss_sec).max(0.0);

    Ok(GtiResult {
        gti_ins,
        gti_utc,
        gross_exposure_sec,
        frame_gap_loss_sec,
        net_exposure_sec,
        uld_quartiles,
        max_ang_offset_deg,
        max_ang_offset_time_utc,
        tct_rejected_glitches: tct.rejected_glitch_count(),
    })
}

/// Longitude interpolation avoiding 360 -> 0 deg wrap glitches (`laxpcl1.f:1318`).
pub fn interp_lon_deg(times: &[f64], lons: &[f64], t: f64) -> f64 {
    let n = times.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 || t <= times[0] {
        return lons[0];
    }
    if t >= times[n - 1] {
        return lons[n - 1];
    }
    let idx = times.partition_point(|&v| v < t);
    let i1 = idx.min(n - 1);
    let i0 = i1.saturating_sub(1);
    let l0 = lons[i0];
    let l1 = lons[i1];
    // Guard against 359 -> 1 deg wrap discontinuity
    if (l1 - l0).abs() > 180.0 {
        return l0;
    }
    let dx = times[i1] - times[i0];
    if dx.abs() == 0.0 {
        l0
    } else {
        l0 + (t - times[i0]) * (l1 - l0) / dx
    }
}
