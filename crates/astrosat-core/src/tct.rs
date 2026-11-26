//! Time Correlation Table (`.tct`) readers and interpolators for AstroSat CZTI and LAXPC.

use std::path::Path;

use astrosat_fits::FitsFile;
use rayon::prelude::*;

use crate::error::{CoreError, Result};

/// CZTI 20-bit second counter wrap modulus (`2^20 = 1048576.0 s`).
pub const CZTI_CLOCK_WRAP_SEC: f64 = 1_048_576.0;
/// Default CZTI Front-End Board (FEB) lag subtracted from `cztseccnt`.
pub const CZTI_DEFAULT_FEB_LAG: f64 = 2.0;
/// LAXPC valid `dt_sps / dt_ins` lower bound (exclusive).
pub const LAXPC_TCT_RATIO_MIN: f64 = 1.000;
/// LAXPC valid `dt_sps / dt_ins` upper bound (inclusive).
pub const LAXPC_TCT_RATIO_MAX: f64 = 1.001;

/// Raw AstroSat Time Correlation Table pairing SPS/UTC seconds with instrument clock seconds.
#[derive(Debug, Clone)]
pub struct TctTable {
    pub sps_time: Vec<f64>,
    pub inst_time: Vec<f64>,
    pub sps_obt: Option<Vec<f64>>,
}

impl TctTable {
    pub fn from_pairs(sps_time: Vec<f64>, inst_time: Vec<f64>) -> Result<Self> {
        if sps_time.len() != inst_time.len() || sps_time.len() < 2 {
            return Err(CoreError::InvalidTct(format!(
                "TCT requires >= 2 matched rows, got sps={} inst={}",
                sps_time.len(),
                inst_time.len()
            )));
        }
        Ok(Self {
            sps_time,
            inst_time,
            sps_obt: None,
        })
    }

    /// Load a TCT table from an AstroSat `.tct` FITS file, auto-detecting the BINTABLE HDU
    /// and column names (`SPS_TIME`/`UT` and `INSTRUMENT_TIME`/`CZT_TIME`/`LAXPC_TIME`).
    pub fn from_fits(path: impl AsRef<Path>, extname: Option<&str>) -> Result<Self> {
        let fits = FitsFile::open(path)?;
        let hdu = if let Some(name) = extname {
            fits.hdu_by_name(name)?
        } else {
            // Find first BINTABLE HDU (usually HDU 1).
            fits.iter_hdus()
                .find(|h| *h.kind() == astrosat_fits::HduKind::BinTable)
                .ok_or_else(|| CoreError::InvalidTct("no BINTABLE HDU found in TCT file".into()))?
        };
        let bt = hdu.as_bintable()?;

        let sps_col = find_first_col(
            &bt,
            &["SPS_TIME", "UT", "UT_TIME", "UTC_TIME", "TIME", "SPS_UTC"],
        )
        .ok_or_else(|| {
            CoreError::InvalidTct("could not locate SPS/UTC time column in TCT".into())
        })?;

        let inst_col = find_first_col(
            &bt,
            &[
                "INSTRUMENT_TIME",
                "INST_TIME",
                "CZT_TIME",
                "CZTI_TIME",
                "LAXPC_TIME",
                "LXP_TIME",
                "ONBOARD_TIME",
                "INST_OBT",
            ],
        )
        .ok_or_else(|| {
            CoreError::InvalidTct("could not locate instrument time column in TCT".into())
        })?;

        let sps_time = bt.read_col_f64(&sps_col)?;
        let inst_time = bt.read_col_f64(&inst_col)?;
        let sps_obt = find_first_col(&bt, &["SPS_OBT", "OBT"])
            .and_then(|c| bt.read_col_f64(&c).ok());

        if sps_time.len() < 2 {
            return Err(CoreError::InvalidTct(format!(
                "TCT file has only {} rows (need >= 2)",
                sps_time.len()
            )));
        }

        Ok(Self {
            sps_time,
            inst_time,
            sps_obt,
        })
    }

    pub fn len(&self) -> usize {
        self.sps_time.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sps_time.is_empty()
    }
}

fn find_first_col(bt: &astrosat_fits::BinTableRef<'_>, candidates: &[&str]) -> Option<String> {
    for &cand in candidates {
        if let Some(spec) = bt.column_spec(cand) {
            return Some(spec.name.clone());
        }
    }
    None
}

/// Segment of monotonic CZTI clock entries between counter resets or 20-bit wrap-arounds.
#[derive(Debug, Clone)]
struct CztiSegment {
    /// Unwrapped monotonic instrument time values (`inst_time` mod 2^20 + wrap offset).
    inst_unwrapped: Vec<f64>,
    /// Corresponding SPS/UTC times.
    sps_time: Vec<f64>,
}

/// CZTI Time Correlation Table interpolator.
///
/// Handles the `-2.0 s` FEB lag, 20-bit counter wrap (`1048576.0 s`), and onboard clock resets.
#[derive(Debug, Clone)]
pub struct CztiTct {
    segments: Vec<CztiSegment>,
    reset_indices: Vec<usize>,
}

impl CztiTct {
    pub fn from_fits(path: impl AsRef<Path>) -> Result<Self> {
        let tbl = TctTable::from_fits(path, None)?;
        Self::new(tbl)
    }

    /// Construct a CZTI TCT interpolator from a raw [`TctTable`].
    ///
    /// Automatically detects both normal 20-bit rollovers (`mod 1048576`) and non-rollover
    /// onboard clock resets.
    pub fn new(table: TctTable) -> Result<Self> {
        if table.len() < 2 {
            return Err(CoreError::InvalidTct("CZTI TCT requires >= 2 rows".into()));
        }

        // Filter out NaN / non-increasing SPS rows.
        let mut pairs: Vec<(f64, f64)> = table
            .sps_time
            .into_iter()
            .zip(table.inst_time)
            .filter(|(s, i)| s.is_finite() && i.is_finite())
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        pairs.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-9);

        if pairs.len() < 2 {
            return Err(CoreError::InvalidTct(
                "CZTI TCT has fewer than 2 valid rows after cleaning".into(),
            ));
        }

        let mut segments: Vec<CztiSegment> = Vec::new();
        let mut reset_indices: Vec<usize> = Vec::new();

        let mut cur_sps = vec![pairs[0].0];
        let mut cur_raw = vec![pairs[0].1.rem_euclid(CZTI_CLOCK_WRAP_SEC)];
        let mut cur_unwrapped = vec![pairs[0].1.rem_euclid(CZTI_CLOCK_WRAP_SEC)];
        let mut wrap_offset = 0.0f64;

        for (idx, &(sps, inst_orig)) in pairs.iter().enumerate().skip(1) {
            let prev_sps = *cur_sps.last().unwrap();
            let dt_sps = sps - prev_sps;
            if dt_sps <= 0.0 {
                continue;
            }
            let raw = inst_orig.rem_euclid(CZTI_CLOCK_WRAP_SEC);
            let prev_raw = *cur_raw.last().unwrap();
            let mut dt_inst = raw - prev_raw;

            // Check for normal 20-bit wrap-around: raw dropped by nearly 1048576.
            if dt_inst < -0.5 * CZTI_CLOCK_WRAP_SEC {
                let unwrapped_dt = dt_inst + CZTI_CLOCK_WRAP_SEC;
                if (unwrapped_dt - dt_sps).abs() < 0.05 * dt_sps.max(100.0) {
                    wrap_offset += CZTI_CLOCK_WRAP_SEC;
                    dt_inst = unwrapped_dt;
                }
            }

            // Detect genuine onboard clock reset (discontinuity inconsistent with elapsed SPS time).
            let rate_mismatch = (dt_inst - dt_sps).abs() > (0.05 * dt_sps).max(10.0);
            if rate_mismatch || dt_inst <= 0.0 {
                if cur_sps.len() >= 2 {
                    segments.push(CztiSegment {
                        inst_unwrapped: cur_unwrapped,
                        sps_time: cur_sps,
                    });
                }
                reset_indices.push(idx);
                wrap_offset = 0.0;
                cur_sps = vec![sps];
                cur_raw = vec![raw];
                cur_unwrapped = vec![raw];
            } else {
                cur_sps.push(sps);
                cur_raw.push(raw);
                cur_unwrapped.push(raw + wrap_offset);
            }
        }

        if cur_sps.len() >= 2 {
            segments.push(CztiSegment {
                inst_unwrapped: cur_unwrapped,
                sps_time: cur_sps,
            });
        }

        if segments.is_empty() {
            return Err(CoreError::InvalidTct(
                "CZTI TCT yielded no valid monotonic segments".into(),
            ));
        }

        Ok(Self {
            segments,
            reset_indices,
        })
    }

    pub fn num_segments(&self) -> usize {
        self.segments.len()
    }

    pub fn reset_count(&self) -> usize {
        self.reset_indices.len()
    }

    /// Interpolate a CZTI `cztseccnt` value to AstroSat UTC seconds.
    ///
    /// Subtracts `feb_lag` (normally `2.0` s) and handles `mod 1048576.0` wrap-around.
    /// If `reset_segment` is `Some(k)`, uses segment `k` explicitly (for onboard reset flags);
    /// otherwise selects the best matching segment automatically.
    pub fn interpolate(
        &self,
        cztseccnt: f64,
        feb_lag: f64,
        reset_segment: Option<usize>,
    ) -> Result<f64> {
        if !cztseccnt.is_finite() || !feb_lag.is_finite() {
            return Err(CoreError::TimeOutOfRange(cztseccnt));
        }
        let t_mod = (cztseccnt - feb_lag).rem_euclid(CZTI_CLOCK_WRAP_SEC);

        if let Some(seg_idx) = reset_segment {
            let seg = self.segments.get(seg_idx).ok_or_else(|| {
                CoreError::InvalidTct(format!("reset segment {seg_idx} out of bounds"))
            })?;
            return Ok(interp_czti_segment(seg, t_mod));
        }

        // Find the segment whose unwrapped range contains t_mod (possibly shifted by k * 2^20).
        let mut best_dist = f64::INFINITY;
        let mut best_val = 0.0;
        for seg in &self.segments {
            let lo = seg.inst_unwrapped[0];
            let hi = *seg.inst_unwrapped.last().unwrap();
            let k_min = ((lo - t_mod) / CZTI_CLOCK_WRAP_SEC).floor() as i64;
            let k_max = ((hi - t_mod) / CZTI_CLOCK_WRAP_SEC).ceil() as i64;
            for k in k_min..=k_max {
                let candidate = t_mod + (k as f64) * CZTI_CLOCK_WRAP_SEC;
                let dist = if candidate < lo {
                    lo - candidate
                } else if candidate > hi {
                    candidate - hi
                } else {
                    0.0
                };
                if dist < best_dist {
                    best_dist = dist;
                    best_val = linear_interp(&seg.inst_unwrapped, &seg.sps_time, candidate);
                    if dist == 0.0 {
                        return Ok(best_val);
                    }
                }
            }
        }
        Ok(best_val)
    }

    /// Parallel batch interpolation for arrays of CZTI `cztseccnt` values.
    pub fn interpolate_batch(&self, times: &[f64], feb_lag: f64) -> Vec<f64> {
        if times.len() >= 16_384 {
            times
                .par_iter()
                .map(|&t| self.interpolate(t, feb_lag, None).unwrap_or(f64::NAN))
                .collect()
        } else {
            times
                .iter()
                .map(|&t| self.interpolate(t, feb_lag, None).unwrap_or(f64::NAN))
                .collect()
        }
    }
}

fn interp_czti_segment(seg: &CztiSegment, t_mod: f64) -> f64 {
    let lo = seg.inst_unwrapped[0];
    let hi = *seg.inst_unwrapped.last().unwrap();
    let mid = 0.5 * (lo + hi);
    let k = ((mid - t_mod) / CZTI_CLOCK_WRAP_SEC).round();
    let t_unwrapped = t_mod + k * CZTI_CLOCK_WRAP_SEC;
    linear_interp(&seg.inst_unwrapped, &seg.sps_time, t_unwrapped)
}

/// LAXPC Time Correlation Table interpolator.
///
/// Rejects glitches where `dt_sps / dt_ins` falls outside `(1.000, 1.001]` and supports
/// fast O(log N) linear and Neville polynomial interpolation.
#[derive(Debug, Clone)]
pub struct LaxpcTct {
    inst_time: Vec<f64>,
    sps_time: Vec<f64>,
    rejected_glitches: usize,
}

impl LaxpcTct {
    pub fn from_fits(path: impl AsRef<Path>) -> Result<Self> {
        let tbl = TctTable::from_fits(path, None)?;
        Self::new(tbl)
    }

    /// Construct a LAXPC TCT interpolator, rejecting glitch rows where `dt_sps / dt_ins`
    /// is outside `(1.000, 1.001]`.
    pub fn new(table: TctTable) -> Result<Self> {
        if table.len() < 2 {
            return Err(CoreError::InvalidTct("LAXPC TCT requires >= 2 rows".into()));
        }

        let raw: Vec<(f64, f64)> = table
            .inst_time
            .into_iter()
            .zip(table.sps_time)
            .filter(|(i, s)| i.is_finite() && s.is_finite())
            .collect();

        if raw.len() < 2 {
            return Err(CoreError::InvalidTct("LAXPC TCT has < 2 finite rows".into()));
        }

        // Find the first valid anchor pair whose ratio is within (1.000, 1.001].
        let mut start_idx = None;
        for i in 0..(raw.len() - 1) {
            let dt_ins = raw[i + 1].0 - raw[i].0;
            let dt_sps = raw[i + 1].1 - raw[i].1;
            if dt_ins > 0.0 {
                let ratio = dt_sps / dt_ins;
                if ratio > LAXPC_TCT_RATIO_MIN && ratio <= LAXPC_TCT_RATIO_MAX {
                    start_idx = Some(i);
                    break;
                }
            }
        }

        let Some(anchor) = start_idx else {
            return Err(CoreError::InvalidTct(
                "LAXPC TCT has no pair with dt_sps/dt_ins in (1.000, 1.001]".into(),
            ));
        };

        let mut inst_time = Vec::with_capacity(raw.len());
        let mut sps_time = Vec::with_capacity(raw.len());
        let mut rejected_glitches = anchor;

        inst_time.push(raw[anchor].0);
        sps_time.push(raw[anchor].1);

        for &(ins, sps) in &raw[(anchor + 1)..] {
            let prev_ins = *inst_time.last().unwrap();
            let prev_sps = *sps_time.last().unwrap();
            let dt_ins = ins - prev_ins;
            let dt_sps = sps - prev_sps;
            if dt_ins > 0.0 {
                let ratio = dt_sps / dt_ins;
                if ratio > LAXPC_TCT_RATIO_MIN && ratio <= LAXPC_TCT_RATIO_MAX {
                    inst_time.push(ins);
                    sps_time.push(sps);
                    continue;
                }
            }
            rejected_glitches += 1;
        }

        if inst_time.len() < 2 {
            return Err(CoreError::InvalidTct(
                "LAXPC TCT has < 2 valid rows after glitch rejection".into(),
            ));
        }

        Ok(Self {
            inst_time,
            sps_time,
            rejected_glitches,
        })
    }

    pub fn len(&self) -> usize {
        self.inst_time.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inst_time.is_empty()
    }

    pub fn rejected_glitch_count(&self) -> usize {
        self.rejected_glitches
    }

    /// O(log N) binary-search linear interpolation from LAXPC instrument clock to AstroSat UTC.
    #[inline]
    pub fn interpolate(&self, inst_t: f64) -> f64 {
        linear_interp(&self.inst_time, &self.sps_time, inst_t)
    }

    /// O(log N) binary-search Neville polynomial interpolation of order `order` (e.g. `order = 4` for 4-point cubic).
    pub fn interpolate_neville(&self, inst_t: f64, order: usize) -> f64 {
        let n = self.inst_time.len();
        let m = order.clamp(2, n.min(8));
        if m == 2 {
            return self.interpolate(inst_t);
        }

        let pos = self.inst_time.partition_point(|&x| x < inst_t);
        let half = m / 2;
        let start = pos.saturating_sub(half).min(n - m);
        let xs = &self.inst_time[start..start + m];
        let mut p = self.sps_time[start..start + m].to_vec();

        // Shift x by xs[0] for numerical conditioning.
        let x0 = xs[0];
        let x_rel = inst_t - x0;
        let mut xr = [0.0f64; 8];
        for i in 0..m {
            xr[i] = xs[i] - x0;
        }

        for k in 1..m {
            for i in 0..(m - k) {
                let denom = xr[i + k] - xr[i];
                if denom.abs() > 0.0 {
                    p[i] = ((x_rel - xr[i]) * p[i + 1] - (x_rel - xr[i + k]) * p[i]) / denom;
                }
            }
        }
        p[0]
    }

    /// Parallel batch interpolation of instrument timestamps.
    pub fn interpolate_batch(&self, times: &[f64]) -> Vec<f64> {
        if times.len() >= 16_384 {
            times.par_iter().map(|&t| self.interpolate(t)).collect()
        } else {
            times.iter().map(|&t| self.interpolate(t)).collect()
        }
    }
}

/// Fast O(log N) piecewise linear interpolation with linear end-extrapolation.
#[inline]
pub(crate) fn linear_interp(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    debug_assert!(n >= 2 && xs.len() == ys.len());
    let idx = xs.partition_point(|&v| v <= x);
    let i = if idx == 0 {
        0
    } else if idx >= n {
        n - 2
    } else {
        idx - 1
    };
    let dx = xs[i + 1] - xs[i];
    if dx == 0.0 {
        return ys[i];
    }
    ys[i] + (x - xs[i]) * (ys[i + 1] - ys[i]) / dx
}
