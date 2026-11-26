//! Stage 2a: `cztgtigen` — MKF threshold parsing and GTI/BTI intersection.

use serde::{Deserialize, Serialize};

use crate::error::{CztiError, Result};
use crate::types::NUM_QUADS;

/// Single parameter threshold rule from `mkfThresholds.txt`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MkfThresholdRule {
    pub column: String,
    pub enabled: bool,
    pub min: Option<f64>,
    pub max: Option<f64>,
}

/// Configuration for `cztgtigen`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GtiGenConfig {
    /// Minimum duration (seconds) for a surviving GTI interval (default: `1.05 s`).
    pub min_exposure: f64,
    /// MKF parameter threshold rules (defaults to standard CZTI `mkfThresholds.txt`).
    pub thresholds: Vec<MkfThresholdRule>,
}

impl Default for GtiGenConfig {
    fn default() -> Self {
        Self {
            min_exposure: 1.05,
            thresholds: default_mkf_thresholds(),
        }
    }
}

/// Standard CZTI Pipeline v3.0.1 MKF thresholds (`config/mkfThresholds.txt`).
pub fn default_mkf_thresholds() -> Vec<MkfThresholdRule> {
    vec![
        MkfThresholdRule {
            column: "ELV".into(),
            enabled: true,
            min: Some(5.0),
            max: None,
        },
        MkfThresholdRule {
            column: "ANG_OFFSET".into(),
            enabled: true,
            min: None,
            max: Some(0.03),
        },
        MkfThresholdRule {
            column: "Q1_MODEID".into(),
            enabled: true,
            min: Some(0.0),
            max: Some(0.0),
        },
        MkfThresholdRule {
            column: "Q1_VetoCounter".into(),
            enabled: true,
            min: None,
            max: Some(450.0),
        },
        MkfThresholdRule {
            column: "Q2_MODEID".into(),
            enabled: true,
            min: Some(0.0),
            max: Some(0.0),
        },
        MkfThresholdRule {
            column: "Q2_VetoCounter".into(),
            enabled: true,
            min: None,
            max: Some(500.0),
        },
        MkfThresholdRule {
            column: "Q3_MODEID".into(),
            enabled: true,
            min: Some(0.0),
            max: Some(0.0),
        },
        MkfThresholdRule {
            column: "Q3_VetoCounter".into(),
            enabled: true,
            min: None,
            max: Some(500.0),
        },
        MkfThresholdRule {
            column: "Q4_MODEID".into(),
            enabled: true,
            min: Some(0.0),
            max: Some(0.0),
        },
        MkfThresholdRule {
            column: "Q4_VetoCounter".into(),
            enabled: true,
            min: None,
            max: Some(480.0),
        },
    ]
}

/// Parse a CZTI `mkfThresholds.txt` file (format: `PARAM  FLAG(0|1)  MIN|-  MAX|-`).
pub fn parse_mkf_threshold_text(text: &str) -> Vec<MkfThresholdRule> {
    let mut rules = Vec::new();
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let enabled = parts[1] == "1";
        let min = if parts[2] == "-" {
            None
        } else {
            parts[2].parse::<f64>().ok()
        };
        let max = if parts[3] == "-" {
            None
        } else {
            parts[3].parse::<f64>().ok()
        };
        if enabled && (min.is_some() || max.is_some()) {
            rules.push(MkfThresholdRule {
                column: parts[0].to_string(),
                enabled,
                min,
                max,
            });
        }
    }
    rules
}

/// Parse a two-column ASCII user GTI text file (`tstart tstop` per line).
pub fn parse_user_gti_text(text: &str) -> Result<Vec<(f64, f64)>> {
    let mut out = Vec::new();
    for (lineno, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            return Err(CztiError::InvalidConfig(format!(
                "Invalid user GTI line {}: {}",
                lineno + 1,
                line
            )));
        }
        let a: f64 = parts[0].parse().map_err(|_| {
            CztiError::InvalidConfig(format!("Invalid GTI start on line {}", lineno + 1))
        })?;
        let b: f64 = parts[1].parse().map_err(|_| {
            CztiError::InvalidConfig(format!("Invalid GTI stop on line {}", lineno + 1))
        })?;
        if b > a {
            out.push((a, b));
        }
    }
    Ok(normalize_intervals(&out))
}

/// Sort and merge overlapping/touching `(start, stop)` intervals.
pub fn normalize_intervals(raw: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut v: Vec<(f64, f64)> = raw.iter().copied().filter(|&(a, b)| b > a).collect();
    if v.is_empty() {
        return v;
    }
    v.sort_by(|x, y| x.0.total_cmp(&y.0).then_with(|| x.1.total_cmp(&y.1)));
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(v.len());
    let mut cur = v[0];
    for &(a, b) in &v[1..] {
        if a <= cur.1 {
            if b > cur.1 {
                cur.1 = b;
            }
        } else {
            out.push(cur);
            cur = (a, b);
        }
    }
    out.push(cur);
    out
}

/// Intersect two sorted/normalized interval lists.
pub fn intersect_intervals(a: &[(f64, f64)], b: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let a = normalize_intervals(a);
    let b = normalize_intervals(b);
    let mut out = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        let lo = a[i].0.max(b[j].0);
        let hi = a[i].1.min(b[j].1);
        if hi > lo {
            out.push((lo, hi));
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    out
}

/// Subtract bad time intervals `bti` from `gti`.
pub fn subtract_intervals(gti: &[(f64, f64)], bti: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let gti = normalize_intervals(gti);
    let bti = normalize_intervals(bti);
    if bti.is_empty() {
        return gti;
    }
    let mut out = Vec::new();
    let mut j = 0usize;
    for &(g0, g1) in &gti {
        let mut cur = g0;
        while j < bti.len() && bti[j].1 <= cur {
            j += 1;
        }
        let mut k = j;
        while k < bti.len() && bti[k].0 < g1 {
            if bti[k].0 > cur {
                out.push((cur, bti[k].0.min(g1)));
            }
            cur = cur.max(bti[k].1);
            if cur >= g1 {
                break;
            }
            k += 1;
        }
        if cur < g1 {
            out.push((cur, g1));
        }
    }
    out
}

/// Drop intervals strictly shorter than `min_exposure` seconds.
pub fn drop_short_intervals(intervals: &[(f64, f64)], min_exposure: f64) -> Vec<(f64, f64)> {
    intervals
        .iter()
        .copied()
        .filter(|&(a, b)| (b - a) >= min_exposure)
        .collect()
}

/// Convert boolean sample validity flags on an MKF `time` grid into contiguous GTI intervals.
///
/// Mirrors `Mkf::get_filtered_gti`: consecutive good MKF samples separated by `<= max_step` (2.5 s)
/// form a single GTI `[t_start, t_stop]` (single isolated good sample spans `[t - 0.5, t + 0.5]`).
pub fn intervals_from_mkf_mask(times: &[f64], good: &[bool]) -> Vec<(f64, f64)> {
    debug_assert_eq!(times.len(), good.len());
    let n = times.len();
    let mut out = Vec::new();
    let mut start: Option<f64> = None;
    let mut prev_t = 0.0f64;

    for i in 0..n {
        if good[i] {
            match start {
                None => {
                    start = Some(times[i]);
                    prev_t = times[i];
                }
                Some(s) => {
                    if times[i] - prev_t > 2.5 {
                        if prev_t > s {
                            out.push((s, prev_t));
                        }
                        start = Some(times[i]);
                    }
                    prev_t = times[i];
                }
            }
        } else if let Some(s) = start.take() {
            if prev_t > s {
                out.push((s, prev_t));
            }
        }
    }
    if let Some(s) = start {
        if prev_t > s {
            out.push((s, prev_t));
        }
    }
    normalize_intervals(&out)
}

/// Generated GTIs for `GTI` (common) and `Q0_GTI..Q3_GTI`.
#[derive(Debug, Clone, Default)]
pub struct CztiGtiOutput {
    pub common: Vec<(f64, f64)>,
    pub quad: [Vec<(f64, f64)>; NUM_QUADS],
}

/// Compute common and per-quadrant GTIs by intersecting `mkf_gti`, input M0 GTIs, optional BTI
/// exclusion, and optional user GTI, then dropping intervals shorter than `min_exposure`.
pub fn generate_czti_gtis(
    mkf_gti: &[(f64, f64)],
    m0_common_gti: &[(f64, f64)],
    m0_quad_gtis: &[Vec<(f64, f64)>; NUM_QUADS],
    m0_bti: Option<&[(f64, f64)]>,
    user_gti: Option<&[(f64, f64)]>,
    min_exposure: f64,
) -> Result<CztiGtiOutput> {
    if mkf_gti.is_empty() {
        return Err(CztiError::EmptyGti);
    }

    let refine = |base: &[(f64, f64)]| -> Vec<(f64, f64)> {
        let mut g = if base.is_empty() {
            mkf_gti.to_vec()
        } else {
            intersect_intervals(mkf_gti, base)
        };
        if let Some(bti) = m0_bti {
            g = subtract_intervals(&g, bti);
        }
        if let Some(ug) = user_gti {
            g = intersect_intervals(&g, ug);
        }
        drop_short_intervals(&g, min_exposure)
    };

    let common = refine(m0_common_gti);
    let quad = [
        refine(&m0_quad_gtis[0]),
        refine(&m0_quad_gtis[1]),
        refine(&m0_quad_gtis[2]),
        refine(&m0_quad_gtis[3]),
    ];

    Ok(CztiGtiOutput { common, quad })
}
