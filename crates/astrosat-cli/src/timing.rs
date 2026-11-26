//! `astrosat timing <bary_evt.fits>` — pulsar epoch-folding (chi-squared) and Z^2_n
//! periodicity detection utility.

use std::f64::consts::TAU;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use astrosat_core::GtiSet;
use astrosat_fits::{FitsFile, HduKind};
use rayon::prelude::*;
use serde::Serialize;

use crate::output::ascii_sparkline;

#[derive(Debug, Clone)]
pub struct TimingParams {
    pub input: PathBuf,
    pub hdu: Option<String>,
    pub gti_file: Option<PathBuf>,
    pub freq_hz: Option<f64>,
    pub period_sec: Option<f64>,
    pub fdot: f64,
    pub fddot: f64,
    pub epoch: Option<f64>,
    pub freq_min: Option<f64>,
    pub freq_max: Option<f64>,
    pub search_span_hz: f64,
    pub df_hz: Option<f64>,
    pub n_trials: Option<usize>,
    pub harmonics: usize,
    pub nbins: usize,
    pub pi_min: Option<i32>,
    pub pi_max: Option<i32>,
    pub unit_filter: Option<u8>,
    pub layer_filter: Option<u8>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TimingTrial {
    pub freq_hz: f64,
    pub period_sec: f64,
    pub z2_n: f64,
    pub chi2: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct TimingResult {
    pub input: PathBuf,
    pub hdu_used: String,
    pub num_events: usize,
    pub t_start: f64,
    pub t_stop: f64,
    pub t_span_sec: f64,
    pub epoch: f64,
    pub fdot: f64,
    pub fddot: f64,
    pub harmonics: usize,
    pub nbins: usize,
    pub n_trials: usize,
    pub best_freq_hz: f64,
    pub best_period_sec: f64,
    pub max_z2_n: f64,
    pub max_chi2: f64,
    pub reduced_chi2: f64,
    pub single_trial_fap: f64,
    pub trial_fap: f64,
    pub detection_sigma: f64,
    pub pulsed_fraction_rms: f64,
    pub profile_counts: Vec<u64>,
    pub profile_normalized: Vec<f64>,
    pub top_trials: Vec<TimingTrial>,
}

pub fn run_timing(params: &TimingParams) -> Result<TimingResult> {
    if params.harmonics == 0 || params.harmonics > 1024 {
        bail!("--harmonics must be in 1..=1024, got {}", params.harmonics);
    }
    if params.nbins < 2 || params.nbins > 65_536 {
        bail!("--nbins must be in 2..=65536, got {}", params.nbins);
    }
    if !params.fdot.is_finite() {
        bail!("--fdot must be finite, got {}", params.fdot);
    }
    if !params.fddot.is_finite() {
        bail!("--fddot must be finite, got {}", params.fddot);
    }
    if !params.search_span_hz.is_finite() || params.search_span_hz < 0.0 {
        bail!(
            "--search-span must be finite and >= 0, got {}",
            params.search_span_hz
        );
    }
    if let Some(ep) = params.epoch {
        if !ep.is_finite() {
            bail!("--epoch must be finite, got {ep}");
        }
    }
    if let Some(df) = params.df_hz {
        if !df.is_finite() || df <= 0.0 {
            bail!("--df must be finite and > 0, got {df}");
        }
    }

    let (hdu_used, mut times) = load_filtered_event_times(params)?;
    if times.len() < 4 {
        bail!(
            "too few events ({}) in {} (HDU {}) for periodicity analysis",
            times.len(),
            params.input.display(),
            hdu_used
        );
    }

    if let Some(ref gti_path) = params.gti_file {
        let gti = GtiSet::from_fits(gti_path, None)
            .with_context(|| format!("failed to load GTI from {}", gti_path.display()))?;
        times.retain(|&t| gti.contains(t));
        if times.len() < 4 {
            bail!("fewer than 4 events remain after GTI filtering");
        }
    }

    let t_start = times.iter().copied().fold(f64::INFINITY, f64::min);
    let t_stop = times.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let t_span_sec = (t_stop - t_start).max(1e-9);
    let epoch = params.epoch.unwrap_or(0.5 * (t_start + t_stop));

    let freqs = build_frequency_grid(params, t_span_sec)?;
    let n_trials = freqs.len();

    let trials: Vec<TimingTrial> = freqs
        .par_iter()
        .map(|&f| evaluate_trial(&times, epoch, f, params.fdot, params.fddot, params.harmonics, params.nbins))
        .collect();

    let best = trials
        .iter()
        .max_by(|a, b| a.z2_n.total_cmp(&b.z2_n))
        .cloned()
        .context("empty trial grid")?;

    let profile_counts = fold_profile(
        &times,
        epoch,
        best.freq_hz,
        params.fdot,
        params.fddot,
        params.nbins,
    );
    let mean_c = (times.len() as f64) / (params.nbins as f64);
    let profile_normalized: Vec<f64> = profile_counts
        .iter()
        .map(|&c| if mean_c > 0.0 { (c as f64) / mean_c } else { 0.0 })
        .collect();

    let var_norm = profile_normalized
        .iter()
        .map(|&r| (r - 1.0) * (r - 1.0))
        .sum::<f64>()
        / (params.nbins as f64);
    let pulsed_fraction_rms = var_norm.sqrt();

    let single_trial_fap = z2_fap(best.z2_n, params.harmonics);
    let trial_fap = if single_trial_fap < 1e-12 {
        (single_trial_fap * (n_trials as f64)).min(1.0)
    } else {
        1.0 - (1.0 - single_trial_fap).powf(n_trials as f64)
    };
    let detection_sigma = fap_to_sigma(trial_fap.max(1e-300));

    let mut sorted = trials;
    sorted.sort_by(|a, b| b.z2_n.total_cmp(&a.z2_n));
    sorted.truncate(5);

    Ok(TimingResult {
        input: params.input.clone(),
        hdu_used,
        num_events: times.len(),
        t_start,
        t_stop,
        t_span_sec,
        epoch,
        fdot: params.fdot,
        fddot: params.fddot,
        harmonics: params.harmonics,
        nbins: params.nbins,
        n_trials,
        best_freq_hz: best.freq_hz,
        best_period_sec: best.period_sec,
        max_z2_n: best.z2_n,
        max_chi2: best.chi2,
        reduced_chi2: best.chi2 / ((params.nbins - 1) as f64),
        single_trial_fap,
        trial_fap,
        detection_sigma,
        pulsed_fraction_rms,
        profile_counts,
        profile_normalized,
        top_trials: sorted,
    })
}

fn build_frequency_grid(params: &TimingParams, t_span_sec: f64) -> Result<Vec<f64>> {
    const MAX_TRIALS: usize = 2_000_000;
    let nyquist_step = 1.0 / (5.0 * t_span_sec);
    if let (Some(fmin), Some(fmax)) = (params.freq_min, params.freq_max) {
        if !fmin.is_finite() || !fmax.is_finite() || fmax <= fmin || fmin <= 0.0 {
            bail!("--freq-max ({fmax}) must be finite and > --freq-min ({fmin}) > 0");
        }
        let n = if let Some(nt) = params.n_trials {
            if nt == 0 || nt > MAX_TRIALS {
                bail!("--n-trials ({nt}) must be in 1..={MAX_TRIALS}");
            }
            nt
        } else {
            let df = params.df_hz.unwrap_or(nyquist_step).max(1e-12);
            let raw = ((fmax - fmin) / df).ceil();
            if !raw.is_finite() || raw >= MAX_TRIALS as f64 {
                MAX_TRIALS
            } else {
                (raw as usize).clamp(1, MAX_TRIALS)
            }
        };
        if n == 1 {
            return Ok(vec![0.5 * (fmin + fmax)]);
        }
        let step = (fmax - fmin) / ((n - 1) as f64);
        return Ok((0..n).map(|i| fmin + (i as f64) * step).collect());
    }

    let center_f = match (params.freq_hz, params.period_sec) {
        (Some(f), _) => f,
        (None, Some(p)) => {
            if !p.is_finite() || p <= 0.0 {
                bail!("--period must be finite and > 0, got {p}");
            }
            1.0 / p
        }
        (None, None) => {
            bail!("Provide --freq <Hz>, --period <sec>, or --freq-min/--freq-max for timing search.");
        }
    };
    if !center_f.is_finite() || center_f <= 0.0 {
        bail!("frequency must be finite and > 0 Hz, got {center_f}");
    }

    if let Some(nt) = params.n_trials {
        if nt == 0 || nt > MAX_TRIALS {
            bail!("--n-trials ({nt}) must be in 1..={MAX_TRIALS}");
        }
    }

    let span = params.search_span_hz;
    if span == 0.0 && params.n_trials.unwrap_or(1) <= 1 {
        return Ok(vec![center_f]);
    }
    let effective_span = if span > 0.0 { span } else { 5.0 * nyquist_step };
    let fmin = (center_f - effective_span).max(1e-9);
    let fmax = center_f + effective_span;
    if !fmin.is_finite() || !fmax.is_finite() || fmax <= fmin {
        bail!("invalid search frequency range [{fmin}, {fmax}]");
    }
    let n = if let Some(nt) = params.n_trials {
        nt
    } else {
        let df = params.df_hz.unwrap_or(nyquist_step).max(1e-12);
        let raw = ((fmax - fmin) / df).ceil();
        if !raw.is_finite() || raw >= MAX_TRIALS as f64 {
            MAX_TRIALS
        } else {
            (raw as usize).clamp(3, MAX_TRIALS)
        }
    };
    if n == 1 {
        return Ok(vec![center_f]);
    }
    let step = (fmax - fmin) / ((n - 1) as f64);
    Ok((0..n).map(|i| fmin + (i as f64) * step).collect())
}

fn load_filtered_event_times(params: &TimingParams) -> Result<(String, Vec<f64>)> {
    let fits = FitsFile::open(&params.input)
        .with_context(|| format!("failed to open event FITS {}", params.input.display()))?;

    let target_hdus: Vec<(usize, String)> = if let Some(ref name) = params.hdu {
        let h = fits
            .hdu_by_name(name)
            .with_context(|| format!("HDU '{name}' not found in {}", params.input.display()))?;
        vec![(h.index(), name.clone())]
    } else if let Ok(h) = fits.hdu_by_name("EVENTS") {
        vec![(h.index(), "EVENTS".to_string())]
    } else if let Ok(h) = fits.hdu_by_name("LAXPC_EVENTS") {
        vec![(h.index(), "LAXPC_EVENTS".to_string())]
    } else if let Ok(h) = fits.hdu_by_name("LAXPC_MODE_EA") {
        vec![(h.index(), "LAXPC_MODE_EA".to_string())]
    } else if fits.hdu_by_name("Q0").is_ok() {
        ["Q0", "Q1", "Q2", "Q3"]
            .into_iter()
            .filter_map(|q| fits.hdu_by_name(q).ok().map(|h| (h.index(), q.to_string())))
            .collect()
    } else {
        let first_bt = fits
            .iter_hdus()
            .find(|h| {
                *h.kind() == HduKind::BinTable
                    && h.as_bintable()
                        .map(|bt| bt.column_index("TIME").is_some())
                        .unwrap_or(false)
            })
            .map(|h| {
                let label = h
                    .extname()
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| format!("HDU_{}", h.index()));
                (h.index(), label)
            })
            .context("no BINTABLE HDU with a TIME column found in FITS file")?;
        vec![first_bt]
    };

    let mut all_times = Vec::new();
    let mut hdu_labels = Vec::with_capacity(target_hdus.len());
    for (hdu_idx, label) in &target_hdus {
        let hdu = fits.hdu(*hdu_idx)?;
        let bt = hdu.as_bintable()?;
        let t = bt.read_col_f64("TIME")?;

        let chan = if params.pi_min.is_some() || params.pi_max.is_some() {
            Some(
                bt.read_col_i32("PI")
                    .or_else(|_| bt.read_col_i32("CHANNEL"))
                    .or_else(|_| bt.read_col_i32("PHA"))
                    .with_context(|| {
                        format!(
                            "--pi-min/--pi-max requested, but HDU '{label}' has no PI/CHANNEL/PHA column"
                        )
                    })?,
            )
        } else {
            None
        };
        let unit = if params.unit_filter.is_some() {
            Some(
                bt.read_col_u8("LAXPC_No")
                    .or_else(|_| bt.read_col_u8("DETECTOR_ID"))
                    .with_context(|| {
                        format!(
                            "--unit requested, but HDU '{label}' has no LAXPC_No/DETECTOR_ID column"
                        )
                    })?,
            )
        } else {
            None
        };
        let layer = if params.layer_filter.is_some() {
            Some(bt.read_col_u8("Layer").with_context(|| {
                format!("--layer requested, but HDU '{label}' has no Layer column")
            })?)
        } else {
            None
        };

        for idx in 0..t.len() {
            let tv = t[idx];
            if !tv.is_finite() {
                continue;
            }
            if let (Some(ch_vec), Some(lo)) = (&chan, params.pi_min) {
                if ch_vec[idx] < lo {
                    continue;
                }
            }
            if let (Some(ch_vec), Some(hi)) = (&chan, params.pi_max) {
                if ch_vec[idx] > hi {
                    continue;
                }
            }
            if let (Some(u_vec), Some(u_req)) = (&unit, params.unit_filter) {
                if u_vec[idx] != u_req {
                    continue;
                }
            }
            if let (Some(l_vec), Some(l_req)) = (&layer, params.layer_filter) {
                if l_vec[idx] != l_req {
                    continue;
                }
            }
            all_times.push(tv);
        }
        hdu_labels.push(label.clone());
    }

    Ok((hdu_labels.join("+"), all_times))
}

#[inline]
fn pulse_phase(t: f64, epoch: f64, f: f64, fdot: f64, fddot: f64) -> f64 {
    let dt = t - epoch;
    let cycles = f * dt + 0.5 * fdot * dt * dt + (1.0 / 6.0) * fddot * dt * dt * dt;
    cycles.rem_euclid(1.0)
}

fn evaluate_trial(
    times: &[f64],
    epoch: f64,
    freq: f64,
    fdot: f64,
    fddot: f64,
    harmonics: usize,
    nbins: usize,
) -> TimingTrial {
    let mut cos_sum = vec![0.0f64; harmonics];
    let mut sin_sum = vec![0.0f64; harmonics];
    let mut bins = vec![0u64; nbins];
    let nbins_f = nbins as f64;

    for &t in times {
        let ph = pulse_phase(t, epoch, freq, fdot, fddot);
        let b = ((ph * nbins_f) as usize).min(nbins - 1);
        bins[b] += 1;

        let ang = TAU * ph;
        // Chebyshev recurrence for cos(k*ang), sin(k*ang) avoids repeated trig calls
        let (s1, c1) = ang.sin_cos();
        let mut ck = c1;
        let mut sk = s1;
        cos_sum[0] += ck;
        sin_sum[0] += sk;
        for k in 1..harmonics {
            let next_c = ck * c1 - sk * s1;
            let next_s = sk * c1 + ck * s1;
            ck = next_c;
            sk = next_s;
            cos_sum[k] += ck;
            sin_sum[k] += sk;
        }
    }

    let n = times.len() as f64;
    let z2_n = (2.0 / n)
        * cos_sum
            .iter()
            .zip(sin_sum.iter())
            .map(|(&c, &s)| c * c + s * s)
            .sum::<f64>();

    let expected = n / nbins_f;
    let chi2 = bins
        .iter()
        .map(|&c| {
            let d = (c as f64) - expected;
            d * d / expected
        })
        .sum::<f64>();

    TimingTrial {
        freq_hz: freq,
        period_sec: 1.0 / freq,
        z2_n,
        chi2,
    }
}

fn fold_profile(
    times: &[f64],
    epoch: f64,
    freq: f64,
    fdot: f64,
    fddot: f64,
    nbins: usize,
) -> Vec<u64> {
    let mut bins = vec![0u64; nbins];
    let nbins_f = nbins as f64;
    for &t in times {
        let ph = pulse_phase(t, epoch, freq, fdot, fddot);
        let b = ((ph * nbins_f) as usize).min(nbins - 1);
        bins[b] += 1;
    }
    bins
}

/// Exact survival probability P(Z^2_n >= z) = exp(-z/2) * sum_{j=0}^{n-1} (z/2)^j / j!.
///
/// Computed term-by-term in log space so large `z2` and high `harmonics` never trigger `0.0 * Inf = NaN`.
pub fn z2_fap(z2: f64, harmonics: usize) -> f64 {
    if !z2.is_finite() || z2 <= 0.0 {
        return 1.0;
    }
    if harmonics == 0 {
        return 1.0;
    }
    let x = 0.5 * z2;
    let ln_x = x.ln();
    let mut log_term = -x;
    let mut sum = log_term.exp();
    for j in 1..harmonics {
        log_term += ln_x - (j as f64).ln();
        sum += log_term.exp();
    }
    sum.clamp(0.0, 1.0)
}

/// Convert a tail probability `p` in `(0, 1)` to equivalent one-sided Gaussian sigma.
pub fn fap_to_sigma(p: f64) -> f64 {
    if !p.is_finite() || p >= 0.5 {
        return 0.0;
    }
    let p = p.max(1e-300);
    // Acklam's rational approximation for p in [1e-6, 0.5) (< 1.15e-9 error)
    if p >= 0.02425 {
        const A: [f64; 6] = [
            -3.969683028665376e+01,
            2.209460984245205e+02,
            -2.759285104469687e+02,
            1.383_577_518_672_69e2,
            -3.066479806614716e+01,
            2.506628277459239e+00,
        ];
        const B: [f64; 5] = [
            -5.447609879822406e+01,
            1.615858368580409e+02,
            -1.556989798598866e+02,
            6.680131188771972e+01,
            -1.328068155288572e+01,
        ];
        let q = 0.5 - p;
        let r = q * q;
        let num = (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q;
        let den = ((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0;
        return (num / den).max(0.0);
    }
    if p >= 1e-6 {
        const C: [f64; 6] = [
            -7.784894002430293e-03,
            -3.223964580411365e-01,
            -2.400758277161838e+00,
            -2.549732539343734e+00,
            4.374664141464968e+00,
            2.938163982698783e+00,
        ];
        const D: [f64; 4] = [
            7.784695709041462e-03,
            3.224671290700398e-01,
            2.445134137142996e+00,
            3.754408661907416e+00,
        ];
        let q = (-2.0 * p.ln()).sqrt();
        let num = ((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5];
        let den = (((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0;
        return (-num / den).max(0.0);
    }
    // Deep tail (p < 1e-6, z > 4.75): asymptotic + Mills-ratio Newton refinement
    let lp = -p.ln();
    let mut z = (2.0 * lp - (2.0 * std::f64::consts::PI * lp).ln().max(0.0))
        .max(1.0)
        .sqrt();
    for _ in 0..3 {
        let pdf = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
        if pdf <= 0.0 {
            break;
        }
        let tail = pdf / (z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0))));
        if tail <= 0.0 {
            break;
        }
        z += (tail - p) / pdf;
    }
    z.max(0.0)
}

pub fn print_timing_human(res: &TimingResult) {
    println!(
        "Timing Analysis: {} (HDU: {}, {} events, span={:.2} s)",
        res.input.display(),
        res.hdu_used,
        res.num_events,
        res.t_span_sec
    );
    println!(
        "  Trials / Grid      : {} trial(s), Z^2_{} & {}-bin epoch folding (epoch={:.3})",
        res.n_trials, res.harmonics, res.nbins, res.epoch
    );
    println!(
        "  Best Frequency     : {:.9} Hz (Period = {:.9} s)",
        res.best_freq_hz, res.best_period_sec
    );
    println!(
        "  Z^2_{:<2} / Chi^2     : Z^2 = {:.2},  Chi^2 = {:.2} (dof = {}, red_chi2 = {:.2})",
        res.harmonics,
        res.max_z2_n,
        res.max_chi2,
        res.nbins - 1,
        res.reduced_chi2
    );
    println!(
        "  Significance       : FAP = {:.3e} ({:.2} sigma), RMS pulsed frac = {:.4}",
        res.trial_fap, res.detection_sigma, res.pulsed_fraction_rms
    );
    println!(
        "  Folded Profile     : [{}]",
        ascii_sparkline(&res.profile_normalized)
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn z2_fap_never_produces_nan_on_large_z2_and_high_harmonics() {
        for (z2, n) in [(55_000.0, 2), (55_000.0, 500), (1e9, 1000), (1e-12, 1)] {
            let fap = z2_fap(z2, n);
            assert!(fap.is_finite() && (0.0..=1.0).contains(&fap), "z2={z2}, n={n}, fap={fap}");
        }
        assert!(fap_to_sigma(0.0).is_finite());
        assert!((fap_to_sigma(0.4) - 0.253347).abs() < 1e-4);
        assert!((fap_to_sigma(0.15865525) - 1.0).abs() < 1e-4);
        assert!((fap_to_sigma(0.02275013) - 2.0).abs() < 1e-4);
        assert!((fap_to_sigma(0.00134990) - 3.0).abs() < 1e-4);
        assert!((fap_to_sigma(1e-12) - 7.03448).abs() < 1e-3);
        assert!(fap_to_sigma(1e-300) > 35.0);
    }
}
