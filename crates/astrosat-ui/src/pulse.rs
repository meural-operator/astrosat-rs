//! High-speed `rayon`-parallelized pulsar periodicity search (`POST /api/pulse/search`)
//! and multi-band / cross-instrument phase folding & pulse alignment (`POST /api/pulse/multiband`).

use std::f64::consts::TAU;
use std::path::{Path, PathBuf};
use std::time::Instant;

use astrosat_core::GtiSet;
use astrosat_fits::{FitsFile, HduKind};
use astrosat_laxpc::background::ChannelMap;
use astrosat_laxpc::types::{utc_to_mission_day, LaxpcUnit};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::http_util::{resolve_existing_file, ApiError, ApiResult};
use crate::science::{
    detect_laxpc_unit, extract_gtis, find_laxpc_channel_map, header_f64_any, header_str_any,
    open_fits_clean,
};

const MAX_SEARCH_TRIALS: usize = 200_000;
const FOLD_CHUNK_SIZE: usize = 65_536;

// ---------------------------------------------------------------------------
// Pulse Search API
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
pub struct PulseSearchRequest {
    pub path: String,
    #[serde(default)]
    pub hdu: Option<String>,
    #[serde(default)]
    pub freq_hz: Option<f64>,
    #[serde(default)]
    pub period_s: Option<f64>,
    #[serde(default)]
    pub freq_min: Option<f64>,
    #[serde(default)]
    pub freq_max: Option<f64>,
    #[serde(default)]
    pub search_span_hz: Option<f64>,
    #[serde(default)]
    pub fdot: Option<f64>,
    #[serde(default)]
    pub fddot: Option<f64>,
    #[serde(default)]
    pub epoch: Option<f64>,
    #[serde(default)]
    pub harmonics: Option<usize>,
    #[serde(default)]
    pub n_trials: Option<usize>,
    #[serde(default)]
    pub nbins: Option<usize>,
    #[serde(default)]
    pub pi_min: Option<i32>,
    #[serde(default)]
    pub pi_max: Option<i32>,
    #[serde(default)]
    pub unit: Option<u8>,
    #[serde(default)]
    pub layer: Option<u8>,
    #[serde(default)]
    pub gti: Option<String>,
    #[serde(default)]
    pub df_hz: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PeriodogramSeries {
    pub freq_hz: Vec<f64>,
    pub z2_n: Vec<f64>,
    pub chi2: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SimpleProfile {
    pub phase_centers: Vec<f64>,
    pub counts: Vec<f64>,
    pub errors: Vec<f64>,
    pub normalized: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PulseSearchResponse {
    pub path: String,
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
    pub elapsed_ms: u64,
    pub best_freq_hz: f64,
    pub best_period_sec: f64,
    pub max_z2_n: f64,
    pub max_chi2: f64,
    pub reduced_chi2: f64,
    pub single_trial_fap: f64,
    pub trial_fap: f64,
    pub detection_sigma: f64,
    pub pulsed_fraction_rms: f64,
    pub periodogram: PeriodogramSeries,
    pub profile: SimpleProfile,
}

/// Handle `POST /api/pulse/search`.
pub fn handle_pulse_search(
    req: &PulseSearchRequest,
    data_dir: &Path,
) -> ApiResult<PulseSearchResponse> {
    let t0_wall = Instant::now();

    let harmonics = req.harmonics.unwrap_or(2);
    if harmonics == 0 || harmonics > 1024 {
        return Err(ApiError::bad_request(format!(
            "harmonics must be in 1..=1024, got {harmonics}"
        )));
    }
    let nbins = req.nbins.unwrap_or(32);
    if !(2..=65_536).contains(&nbins) {
        return Err(ApiError::bad_request(format!(
            "nbins must be in 2..=65536, got {nbins}"
        )));
    }
    let fdot = req.fdot.unwrap_or(0.0);
    if !fdot.is_finite() {
        return Err(ApiError::bad_request("fdot must be finite"));
    }
    let fddot = req.fddot.unwrap_or(0.0);
    if !fddot.is_finite() {
        return Err(ApiError::bad_request("fddot must be finite"));
    }
    let search_span = req.search_span_hz.unwrap_or(0.0);
    if !search_span.is_finite() || search_span < 0.0 {
        return Err(ApiError::bad_request(
            "search_span_hz must be finite and >= 0",
        ));
    }
    if let Some(ep) = req.epoch {
        if !ep.is_finite() {
            return Err(ApiError::bad_request("epoch must be finite"));
        }
    }
    if let Some(df) = req.df_hz {
        if !df.is_finite() || df <= 0.0 {
            return Err(ApiError::bad_request(format!(
                "df_hz must be finite and > 0, got {df}"
            )));
        }
    }

    let path = resolve_existing_file(&req.path, data_dir)?;
    let (hdu_used, mut times) = load_filtered_event_times(
        &path,
        req.hdu.as_deref(),
        req.pi_min,
        req.pi_max,
        req.unit,
        req.layer,
    )?;

    if times.len() < 4 {
        return Err(ApiError::bad_request(format!(
            "too few events ({}) in {} (HDU {}) for periodicity search",
            times.len(),
            path.display(),
            hdu_used
        )));
    }

    if let Some(gti_raw) = req.gti.as_deref().map(str::trim).filter(|s| !s.is_empty()) {
        let gti_path = resolve_existing_file(gti_raw, data_dir)?;
        let gti = load_gti_set(&gti_path)?;
        times.retain(|&t| gti.contains(t));
        if times.len() < 4 {
            return Err(ApiError::bad_request(
                "fewer than 4 events remain after GTI filtering",
            ));
        }
    }

    let t_start = times.iter().copied().fold(f64::INFINITY, f64::min);
    let t_stop = times.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let t_span_sec = (t_stop - t_start).max(1e-9);
    let epoch = req.epoch.unwrap_or(0.5 * (t_start + t_stop));

    let freqs = build_frequency_grid(
        req.freq_hz,
        req.period_s,
        req.freq_min,
        req.freq_max,
        search_span,
        req.df_hz,
        req.n_trials,
        t_span_sec,
    )?;
    let n_trials = freqs.len();

    let trials: Vec<(f64, f64, f64)> = freqs
        .par_iter()
        .map(|&f| {
            let (z2, chi2) = evaluate_trial(&times, epoch, f, fdot, fddot, harmonics, nbins);
            (f, z2, chi2)
        })
        .collect();

    let (best_freq_hz, max_z2_n, max_chi2) = trials
        .iter()
        .copied()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .ok_or_else(|| ApiError::bad_request("empty frequency trial grid"))?;

    let counts_u64 = fold_counts_parallel(&times, epoch, best_freq_hz, fdot, fddot, nbins);
    let counts: Vec<f64> = counts_u64.iter().map(|&c| c as f64).collect();
    let errors: Vec<f64> = counts.iter().map(|&c| c.sqrt()).collect();
    let mean_c = (times.len() as f64) / (nbins as f64);
    let normalized: Vec<f64> = counts
        .iter()
        .map(|&c| if mean_c > 0.0 { c / mean_c } else { 0.0 })
        .collect();
    let pulsed_fraction_rms = (normalized
        .iter()
        .map(|&r| (r - 1.0) * (r - 1.0))
        .sum::<f64>()
        / (nbins as f64))
        .sqrt();

    let single_trial_fap = z2_fap(max_z2_n, harmonics);
    let trial_fap = if single_trial_fap < 1e-12 {
        (single_trial_fap * (n_trials as f64)).min(1.0)
    } else {
        1.0 - (1.0 - single_trial_fap).powf(n_trials as f64)
    };
    let detection_sigma = fap_to_sigma(trial_fap.max(1e-300));

    let phase_centers: Vec<f64> = (0..nbins)
        .map(|i| ((i as f64) + 0.5) / (nbins as f64))
        .collect();

    let periodogram = PeriodogramSeries {
        freq_hz: trials.iter().map(|t| t.0).collect(),
        z2_n: trials.iter().map(|t| t.1).collect(),
        chi2: trials.iter().map(|t| t.2).collect(),
    };

    Ok(PulseSearchResponse {
        path: path.to_string_lossy().into_owned(),
        hdu_used,
        num_events: times.len(),
        t_start,
        t_stop,
        t_span_sec,
        epoch,
        fdot,
        fddot,
        harmonics,
        nbins,
        n_trials,
        elapsed_ms: t0_wall.elapsed().as_millis() as u64,
        best_freq_hz,
        best_period_sec: 1.0 / best_freq_hz,
        max_z2_n,
        max_chi2,
        reduced_chi2: max_chi2 / ((nbins - 1) as f64),
        single_trial_fap,
        trial_fap,
        detection_sigma,
        pulsed_fraction_rms,
        periodogram,
        profile: SimpleProfile {
            phase_centers,
            counts,
            errors,
            normalized,
        },
    })
}

// ---------------------------------------------------------------------------
// Multi-Band Phase Folding & Alignment API
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
pub struct BandSpec {
    pub label: String,
    #[serde(default)]
    pub emin_kev: Option<f64>,
    #[serde(default)]
    pub emax_kev: Option<f64>,
    #[serde(default)]
    pub pi_min: Option<i32>,
    #[serde(default)]
    pub pi_max: Option<i32>,
    #[serde(default)]
    pub chan_min: Option<i32>,
    #[serde(default)]
    pub chan_max: Option<i32>,
    #[serde(default)]
    pub unit: Option<u8>,
    #[serde(default)]
    pub layer: Option<u8>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct MultibandInputSpec {
    #[serde(default)]
    pub label: Option<String>,
    pub path: String,
    #[serde(default)]
    pub hdu: Option<String>,
    #[serde(default)]
    pub caldb: Option<String>,
    #[serde(default)]
    pub use_weights: Option<bool>,
    pub bands: Vec<BandSpec>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct MultibandRequest {
    pub inputs: Vec<MultibandInputSpec>,
    pub freq_hz: f64,
    #[serde(default)]
    pub fdot: Option<f64>,
    #[serde(default)]
    pub fddot: Option<f64>,
    #[serde(default)]
    pub epoch: Option<f64>,
    #[serde(default)]
    pub nbins: Option<usize>,
    #[serde(default)]
    pub harmonics: Option<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct BandProfileResult {
    pub input_label: String,
    pub band_label: String,
    pub full_label: String,
    pub path: String,
    pub hdu_used: String,
    pub weighted: bool,
    pub num_events: usize,
    pub total_weight: f64,
    pub exposure_sec: f64,
    pub counts: Vec<f64>,
    pub errors: Vec<f64>,
    pub rate: Vec<f64>,
    pub rate_errors: Vec<f64>,
    pub normalized: Vec<f64>,
    pub normalized_errors: Vec<f64>,
    pub pulsed_fraction_minmax: f64,
    pub pulsed_fraction_rms: f64,
    pub z2_n: f64,
    pub chi2: f64,
    pub reduced_chi2: f64,
    pub peak_phase: f64,
    pub phase_lag_cycles: f64,
    pub phase_lag_ms: f64,
    pub xcorr_peak: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct MultibandResponse {
    pub freq_hz: f64,
    pub period_sec: f64,
    pub fdot: f64,
    pub fddot: f64,
    pub epoch: f64,
    pub nbins: usize,
    pub harmonics: usize,
    pub elapsed_ms: u64,
    pub reference_band: Option<String>,
    pub phase_centers: Vec<f64>,
    pub bands: Vec<BandProfileResult>,
}

struct LoadedEventInput {
    input_label: String,
    path: PathBuf,
    hdu_used: String,
    weighted: bool,
    exposure_sec: f64,
    times: Vec<f64>,
    energies_kev: Option<Vec<f32>>,
    channels: Option<Vec<i32>>,
    units: Option<Vec<u8>>,
    layers: Option<Vec<u8>>,
    weights: Option<Vec<f64>>,
    bands: Vec<BandSpec>,
}

/// Handle `POST /api/pulse/multiband`.
pub fn handle_pulse_multiband(
    req: &MultibandRequest,
    data_dir: &Path,
) -> ApiResult<MultibandResponse> {
    let t0_wall = Instant::now();

    if !req.freq_hz.is_finite() || req.freq_hz <= 0.0 {
        return Err(ApiError::bad_request(format!(
            "freq_hz must be finite and > 0, got {}",
            req.freq_hz
        )));
    }
    let fdot = req.fdot.unwrap_or(0.0);
    if !fdot.is_finite() {
        return Err(ApiError::bad_request("fdot must be finite"));
    }
    let fddot = req.fddot.unwrap_or(0.0);
    if !fddot.is_finite() {
        return Err(ApiError::bad_request("fddot must be finite"));
    }
    if let Some(ep) = req.epoch {
        if !ep.is_finite() {
            return Err(ApiError::bad_request("epoch must be finite"));
        }
    }
    let nbins = req.nbins.unwrap_or(64);
    if !(2..=65_536).contains(&nbins) {
        return Err(ApiError::bad_request(format!(
            "nbins must be in 2..=65536, got {nbins}"
        )));
    }
    let harmonics = req.harmonics.unwrap_or(2);
    if harmonics == 0 || harmonics > 1024 {
        return Err(ApiError::bad_request(format!(
            "harmonics must be in 1..=1024, got {harmonics}"
        )));
    }
    if req.inputs.is_empty() {
        return Err(ApiError::bad_request("inputs must contain at least one file"));
    }

    // Validate band parameters before loading files
    for inp in &req.inputs {
        if inp.bands.is_empty() {
            return Err(ApiError::bad_request(format!(
                "input '{}' must specify at least one band",
                inp.path
            )));
        }
        for b in &inp.bands {
            if let Some(v) = b.emin_kev {
                if !v.is_finite() {
                    return Err(ApiError::bad_request("emin_kev must be finite"));
                }
            }
            if let Some(v) = b.emax_kev {
                if !v.is_finite() {
                    return Err(ApiError::bad_request("emax_kev must be finite"));
                }
            }
            if let (Some(lo), Some(hi)) = (b.emin_kev, b.emax_kev) {
                if hi < lo {
                    return Err(ApiError::bad_request(format!(
                        "band '{}': emax_kev ({hi}) < emin_kev ({lo})",
                        b.label
                    )));
                }
            }
            if let (Some(lo), Some(hi)) = (b.pi_min, b.pi_max) {
                if hi < lo {
                    return Err(ApiError::bad_request(format!(
                        "band '{}': pi_max ({hi}) < pi_min ({lo})",
                        b.label
                    )));
                }
            }
            if let (Some(lo), Some(hi)) = (b.chan_min, b.chan_max) {
                if hi < lo {
                    return Err(ApiError::bad_request(format!(
                        "band '{}': chan_max ({hi}) < chan_min ({lo})",
                        b.label
                    )));
                }
            }
        }
    }

    // Load all input FITS event arrays
    let mut loaded_inputs = Vec::with_capacity(req.inputs.len());
    let mut global_tmin = f64::INFINITY;
    let mut global_tmax = f64::NEG_INFINITY;

    for inp in &req.inputs {
        let loaded = load_multiband_input(inp, data_dir)?;
        for &t in &loaded.times {
            if t < global_tmin {
                global_tmin = t;
            }
            if t > global_tmax {
                global_tmax = t;
            }
        }
        loaded_inputs.push(loaded);
    }

    let epoch = req.epoch.unwrap_or_else(|| {
        if global_tmin.is_finite() && global_tmax.is_finite() {
            0.5 * (global_tmin + global_tmax)
        } else {
            0.0
        }
    });

    let period_sec = 1.0 / req.freq_hz;
    let phase_centers: Vec<f64> = (0..nbins)
        .map(|i| ((i as f64) + 0.5) / (nbins as f64))
        .collect();

    // Fold every input's bands in parallel using single-pass rayon chunk reduction
    let mut band_results: Vec<BandProfileResult> = Vec::new();
    for loaded in &loaded_inputs {
        let folded = fold_input_bands_parallel(
            loaded,
            epoch,
            req.freq_hz,
            fdot,
            fddot,
            nbins,
            harmonics,
        );
        band_results.extend(folded);
    }

    // Compute cross-correlation phase lag relative to the first non-empty reference band
    let ref_idx = band_results.iter().position(|b| b.num_events > 0);
    let reference_band = ref_idx.map(|i| band_results[i].full_label.clone());

    if let Some(ri) = ref_idx {
        let ref_profile = band_results[ri].counts.clone();
        for (idx, b) in band_results.iter_mut().enumerate() {
            if idx == ri {
                b.phase_lag_cycles = 0.0;
                b.phase_lag_ms = 0.0;
                b.xcorr_peak = 1.0;
            } else if b.num_events > 0 {
                let (lag_cycles, xcorr_peak) =
                    circular_xcorr_phase_lag(&ref_profile, &b.counts);
                b.phase_lag_cycles = lag_cycles;
                b.phase_lag_ms = lag_cycles * period_sec * 1000.0;
                b.xcorr_peak = xcorr_peak;
            }
        }
    }

    Ok(MultibandResponse {
        freq_hz: req.freq_hz,
        period_sec,
        fdot,
        fddot,
        epoch,
        nbins,
        harmonics,
        elapsed_ms: t0_wall.elapsed().as_millis() as u64,
        reference_band,
        phase_centers,
        bands: band_results,
    })
}

fn load_multiband_input(
    spec: &MultibandInputSpec,
    data_dir: &Path,
) -> ApiResult<LoadedEventInput> {
    let path = resolve_existing_file(&spec.path, data_dir)?;
    let fits = open_fits_clean(&path)?;
    let target_hdus = select_event_hdus(&fits, spec.hdu.as_deref(), &path)?;

    let need_energy = spec
        .bands
        .iter()
        .any(|b| b.emin_kev.is_some() || b.emax_kev.is_some());
    let need_chan = spec
        .bands
        .iter()
        .any(|b| b.pi_min.is_some() || b.pi_max.is_some() || b.chan_min.is_some() || b.chan_max.is_some());
    let need_unit = spec.bands.iter().any(|b| b.unit.is_some());
    let need_layer = spec.bands.iter().any(|b| b.layer.is_some());
    let want_weights = spec.use_weights.unwrap_or(false);

    let instrume = header_str_any(&fits, target_hdus[0].0, "INSTRUME");
    let detnam = header_str_any(&fits, target_hdus[0].0, "DETNAM");
    let tstart_hdr = header_f64_any(&fits, target_hdus[0].0, "TSTART");
    let exposure_sec = compute_event_exposure(&fits, &target_hdus);

    let fname_lo = path
        .file_name()
        .map(|s| s.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    let inst_up = instrume.as_deref().unwrap_or("").to_ascii_uppercase();
    let is_laxpc = inst_up.contains("LAXPC")
        || fname_lo.contains("lxp")
        || fname_lo.contains("laxpc")
        || target_hdus
            .iter()
            .any(|(_, name)| name.to_ascii_uppercase().contains("LAXPC"));
    let is_czti = inst_up.contains("CZTI")
        || fname_lo.contains("czt")
        || target_hdus.iter().any(|(_, name)| matches!(name.as_str(), "Q0" | "Q1" | "Q2" | "Q3"));

    let caldb_explicit = spec.caldb.as_deref().map(PathBuf::from);
    let ndays = tstart_hdr.map(utc_to_mission_day).unwrap_or(100.0);
    let laxpc_maps: [Option<ChannelMap>; 3] = if need_energy && is_laxpc {
        [
            find_laxpc_channel_map(LaxpcUnit::Lxp10, ndays, caldb_explicit.as_deref(), data_dir)
                .map(|(m, _)| m),
            find_laxpc_channel_map(LaxpcUnit::Lxp20, ndays, caldb_explicit.as_deref(), data_dir)
                .map(|(m, _)| m),
            find_laxpc_channel_map(LaxpcUnit::Lxp30, ndays, caldb_explicit.as_deref(), data_dir)
                .map(|(m, _)| m),
        ]
    } else {
        [None, None, None]
    };
    let default_laxpc_unit = detect_laxpc_unit(detnam.as_deref(), &fname_lo).unwrap_or(LaxpcUnit::Lxp10);

    let mut all_times = Vec::new();
    let mut all_energies: Option<Vec<f32>> = if need_energy { Some(Vec::new()) } else { None };
    let mut all_channels: Option<Vec<i32>> = if need_chan { Some(Vec::new()) } else { None };
    let mut all_units: Option<Vec<u8>> = if need_unit { Some(Vec::new()) } else { None };
    let mut all_layers: Option<Vec<u8>> = if need_layer { Some(Vec::new()) } else { None };
    let mut all_weights: Option<Vec<f64>> = None;
    let mut any_weights = false;
    let mut labels = Vec::with_capacity(target_hdus.len());

    for (hdu_idx, label) in &target_hdus {
        let hdu = fits
            .hdu(*hdu_idx)
            .map_err(|e| ApiError::bad_request(format!("failed reading HDU '{label}': {e}")))?;
        let hdr = hdu.header();
        let bt = hdu
            .as_bintable()
            .map_err(|e| ApiError::bad_request(format!("HDU '{label}' is not a BINTABLE: {e}")))?;

        let times = bt
            .read_col_f64("TIME")
            .map_err(|e| ApiError::bad_request(format!("HDU '{label}' missing TIME: {e}")))?;
        let n = times.len();

        let chan_col_name = if bt.column_index("PI").is_some() {
            Some("PI")
        } else if bt.column_index("CHANNEL").is_some() {
            Some("CHANNEL")
        } else if bt.column_index("PHA").is_some() {
            Some("PHA")
        } else {
            None
        };

        let chan_vec = if need_chan || (need_energy && bt.column_index("ENERGY").is_none()) {
            let cname = chan_col_name.ok_or_else(|| {
                ApiError::bad_request(format!(
                    "HDU '{label}' has no ENERGY, PI, CHANNEL, or PHA column for band filtering"
                ))
            })?;
            Some(bt.read_col_i32(cname).map_err(|e| {
                ApiError::bad_request(format!("failed reading {cname} in '{label}': {e}"))
            })?)
        } else {
            None
        };

        let unit_vec = if need_unit || (need_energy && is_laxpc) {
            bt.read_col_u8("LAXPC_No")
                .or_else(|_| bt.read_col_u8("DETECTOR_ID"))
                .ok()
        } else {
            None
        };
        if need_unit && unit_vec.is_none() {
            return Err(ApiError::bad_request(format!(
                "unit filter requested, but HDU '{label}' has no LAXPC_No/DETECTOR_ID column"
            )));
        }

        let layer_vec = if need_layer {
            Some(bt.read_col_u8("Layer").map_err(|_| {
                ApiError::bad_request(format!(
                    "layer filter requested, but HDU '{label}' has no Layer column"
                ))
            })?)
        } else {
            None
        };

        let energy_vec: Option<Vec<f32>> = if need_energy {
            if bt.column_index("ENERGY").is_some() {
                let e64 = bt.read_col_f64("ENERGY").map_err(|e| {
                    ApiError::bad_request(format!("failed reading ENERGY in '{label}': {e}"))
                })?;
                Some(e64.into_iter().map(|v| v as f32).collect())
            } else if is_laxpc {
                let chans = chan_vec.as_ref().expect("chan_vec loaded above");
                // Determine channel scale from TLMAXn header on the channel column
                let chan_scale = infer_laxpc_channel_scale(hdr, &bt, chan_col_name.unwrap(), chans);
                let luts: [Option<Vec<f32>>; 3] = [
                    laxpc_maps[0]
                        .as_ref()
                        .map(|m| build_laxpc_energy_lut(m, chan_scale)),
                    laxpc_maps[1]
                        .as_ref()
                        .map(|m| build_laxpc_energy_lut(m, chan_scale)),
                    laxpc_maps[2]
                        .as_ref()
                        .map(|m| build_laxpc_energy_lut(m, chan_scale)),
                ];
                let fallback_lut = luts[default_laxpc_unit.index()]
                    .as_ref()
                    .or_else(|| luts.iter().flatten().next())
                    .ok_or_else(|| {
                        ApiError::bad_request(
                            "LAXPC keV energy band requested, but no lxp{10,20,30}.channel CALDB map was found (pass `caldb` or use `chan_min`/`chan_max`)",
                        )
                    })?;
                let mut e_out = Vec::with_capacity(n);
                for i in 0..n {
                    let u_idx = unit_vec
                        .as_ref()
                        .and_then(|uv| {
                            let id = uv[i];
                            if (1..=3).contains(&id) {
                                Some((id - 1) as usize)
                            } else {
                                None
                            }
                        })
                        .unwrap_or_else(|| default_laxpc_unit.index());
                    let lut = luts[u_idx].as_ref().unwrap_or(fallback_lut);
                    let c = chans[i].max(0) as usize;
                    e_out.push(lut[c.min(lut.len() - 1)]);
                }
                Some(e_out)
            } else if is_czti {
                let chans = chan_vec.as_ref().expect("chan_vec loaded above");
                // Standard CZTI PI mapping: E_mid = 5.25 + 0.5 * PI
                Some(
                    chans
                        .iter()
                        .map(|&pi| (5.25 + 0.5 * (pi as f64)) as f32)
                        .collect(),
                )
            } else {
                return Err(ApiError::bad_request(format!(
                    "HDU '{label}' has no ENERGY column and is neither LAXPC nor CZTI; use pi_min/pi_max or chan_min/chan_max"
                )));
            }
        } else {
            None
        };

        let weight_vec = if want_weights && bt.column_index("WEIGHT").is_some() {
            any_weights = true;
            Some(
                bt.read_col_f64("WEIGHT")
                    .unwrap_or_else(|_| vec![1.0; n]),
            )
        } else {
            None
        };

        // Filter non-finite timestamps consistently across all columns
        for i in 0..n {
            let tv = times[i];
            if !tv.is_finite() {
                continue;
            }
            all_times.push(tv);
            if let (Some(dst), Some(src)) = (&mut all_energies, &energy_vec) {
                dst.push(src[i]);
            }
            if let (Some(dst), Some(src)) = (&mut all_channels, &chan_vec) {
                dst.push(src[i]);
            }
            if let (Some(dst), Some(src)) = (&mut all_units, &unit_vec) {
                dst.push(src[i]);
            }
            if let (Some(dst), Some(src)) = (&mut all_layers, &layer_vec) {
                dst.push(src[i]);
            }
            if any_weights {
                let w_dst = all_weights.get_or_insert_with(|| vec![1.0; all_times.len() - 1]);
                let wv = weight_vec
                    .as_ref()
                    .map(|wv| if wv[i].is_finite() { wv[i] } else { 0.0 })
                    .unwrap_or(1.0);
                w_dst.push(wv);
            }
        }

        labels.push(label.clone());
    }

    let input_label = spec
        .label
        .as_deref()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .or_else(|| instrume.clone())
        .unwrap_or_else(|| {
            path.file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "Input".to_string())
        });

    Ok(LoadedEventInput {
        input_label,
        path,
        hdu_used: labels.join("+"),
        weighted: any_weights,
        exposure_sec,
        times: all_times,
        energies_kev: all_energies,
        channels: all_channels,
        units: all_units,
        layers: all_layers,
        weights: all_weights,
        bands: spec.bands.clone(),
    })
}

fn infer_laxpc_channel_scale(
    hdr: &astrosat_fits::Header,
    bt: &astrosat_fits::BinTableRef<'_>,
    col_name: &str,
    chans: &[i32],
) -> usize {
    if let Some(ci) = bt.column_index(col_name) {
        let tlmax_key = format!("TLMAX{}", ci + 1);
        if let Some(tlmax) = hdr.get_i64(&tlmax_key) {
            if tlmax == 255 {
                return 256;
            } else if tlmax == 511 {
                return 512;
            } else if tlmax == 1023 {
                return 1024;
            } else if tlmax > 0 && tlmax <= 4096 {
                return (tlmax + 1) as usize;
            }
        }
    }
    let max_ch = chans.iter().copied().max().unwrap_or(511);
    if max_ch >= 512 {
        1024
    } else if max_ch >= 256 {
        512
    } else {
        256
    }
}

fn build_laxpc_energy_lut(cmap: &ChannelMap, nch: usize) -> Vec<f32> {
    let map_len = cmap.len().max(1);
    let mut lut = vec![0.0f32; nch.max(1)];
    if nch == map_len {
        for (i, slot) in lut.iter_mut().enumerate() {
            *slot = (0.5 * (cmap.e_min_kev[i] + cmap.e_max_kev[i])) as f32;
        }
    } else if nch > map_len && nch.is_multiple_of(map_len) {
        let factor = nch / map_len;
        let factor_f = factor as f64;
        for (c, slot) in lut.iter_mut().enumerate() {
            let native = (c / factor).min(map_len - 1);
            let sub = (c % factor) as f64;
            let lo0 = cmap.e_min_kev[native];
            let hi0 = cmap.e_max_kev[native];
            let de = (hi0 - lo0) / factor_f;
            *slot = (lo0 + (sub + 0.5) * de) as f32;
        }
    } else if map_len > nch && map_len.is_multiple_of(nch.max(1)) {
        let factor = map_len / nch.max(1);
        for (c, slot) in lut.iter_mut().enumerate() {
            let lo_idx = (c * factor).min(map_len - 1);
            let hi_idx = ((c + 1) * factor - 1).min(map_len - 1);
            *slot = (0.5 * (cmap.e_min_kev[lo_idx] + cmap.e_max_kev[hi_idx])) as f32;
        }
    } else {
        let ratio = (map_len as f64) / (nch.max(1) as f64);
        for (c, slot) in lut.iter_mut().enumerate() {
            let idx = (((c as f64) * ratio).floor() as usize).min(map_len - 1);
            *slot = (0.5 * (cmap.e_min_kev[idx] + cmap.e_max_kev[idx])) as f32;
        }
    }
    lut
}

fn compute_event_exposure(fits: &FitsFile, target_hdus: &[(usize, String)]) -> f64 {
    if let Some(exp) = header_f64_any(fits, target_hdus[0].0, "EXPOSURE") {
        if exp > 0.0 {
            return exp;
        }
    }
    let gtis = extract_gtis(fits);
    if !gtis.is_empty() {
        let gti_exp: f64 = gtis.iter().map(|g| (g.stop - g.start).max(0.0)).sum();
        if gti_exp > 0.0 {
            return gti_exp;
        }
    }
    if let (Some(ts), Some(te)) = (
        header_f64_any(fits, target_hdus[0].0, "TSTART"),
        header_f64_any(fits, target_hdus[0].0, "TSTOP"),
    ) {
        if te > ts {
            return te - ts;
        }
    }
    0.0
}

#[derive(Clone)]
struct BandAccum {
    num_events: usize,
    total_weight: f64,
    total_w2: f64,
    bins: Vec<f64>,
    bins_w2: Vec<f64>,
    cos_sum: Vec<f64>,
    sin_sum: Vec<f64>,
}

impl BandAccum {
    fn new(nbins: usize, harmonics: usize) -> Self {
        Self {
            num_events: 0,
            total_weight: 0.0,
            total_w2: 0.0,
            bins: vec![0.0; nbins],
            bins_w2: vec![0.0; nbins],
            cos_sum: vec![0.0; harmonics],
            sin_sum: vec![0.0; harmonics],
        }
    }

    fn merge(&mut self, other: &Self) {
        self.num_events += other.num_events;
        self.total_weight += other.total_weight;
        self.total_w2 += other.total_w2;
        for (a, &b) in self.bins.iter_mut().zip(other.bins.iter()) {
            *a += b;
        }
        for (a, &b) in self.bins_w2.iter_mut().zip(other.bins_w2.iter()) {
            *a += b;
        }
        for (a, &b) in self.cos_sum.iter_mut().zip(other.cos_sum.iter()) {
            *a += b;
        }
        for (a, &b) in self.sin_sum.iter_mut().zip(other.sin_sum.iter()) {
            *a += b;
        }
    }
}

fn fold_input_bands_parallel(
    loaded: &LoadedEventInput,
    epoch: f64,
    freq_hz: f64,
    fdot: f64,
    fddot: f64,
    nbins: usize,
    harmonics: usize,
) -> Vec<BandProfileResult> {
    let nb = loaded.bands.len();
    let n_ev = loaded.times.len();
    let nbins_f = nbins as f64;

    let n_chunks = n_ev.div_ceil(FOLD_CHUNK_SIZE).max(1);
    let chunk_accums: Vec<Vec<BandAccum>> = (0..n_chunks)
        .into_par_iter()
        .map(|chunk_idx| {
            let start = chunk_idx * FOLD_CHUNK_SIZE;
            let end = (start + FOLD_CHUNK_SIZE).min(n_ev);
            let mut local: Vec<BandAccum> = (0..nb)
                .map(|_| BandAccum::new(nbins, harmonics))
                .collect();

            for i in start..end {
                let t = loaded.times[i];
                let e_opt = loaded.energies_kev.as_ref().map(|v| v[i] as f64);
                let ch_opt = loaded.channels.as_ref().map(|v| v[i]);
                let u_opt = loaded.units.as_ref().map(|v| v[i]);
                let l_opt = loaded.layers.as_ref().map(|v| v[i]);
                let w = loaded.weights.as_ref().map(|v| v[i]).unwrap_or(1.0);

                // Lazily compute phase & trig only if at least one band matches this event
                let mut phase_cache: Option<(usize, f64, f64)> = None;

                for (bi, band) in loaded.bands.iter().enumerate() {
                    if !event_matches_band(band, e_opt, ch_opt, u_opt, l_opt) {
                        continue;
                    }
                    let (bin_idx, c1, s1) = *phase_cache.get_or_insert_with(|| {
                        let ph = pulse_phase(t, epoch, freq_hz, fdot, fddot);
                        let b = ((ph * nbins_f) as usize).min(nbins - 1);
                        let (s1, c1) = (TAU * ph).sin_cos();
                        (b, c1, s1)
                    });

                    let acc = &mut local[bi];
                    acc.num_events += 1;
                    acc.total_weight += w;
                    acc.total_w2 += w * w;
                    acc.bins[bin_idx] += w;
                    acc.bins_w2[bin_idx] += w * w;

                    let mut ck = c1;
                    let mut sk = s1;
                    acc.cos_sum[0] += w * ck;
                    acc.sin_sum[0] += w * sk;
                    for k in 1..harmonics {
                        let next_c = ck * c1 - sk * s1;
                        let next_s = sk * c1 + ck * s1;
                        ck = next_c;
                        sk = next_s;
                        acc.cos_sum[k] += w * ck;
                        acc.sin_sum[k] += w * sk;
                    }
                }
            }
            local
        })
        .collect();

    // Deterministic ordered merge across chunks
    let mut merged: Vec<BandAccum> = (0..nb)
        .map(|_| BandAccum::new(nbins, harmonics))
        .collect();
    for chunk in &chunk_accums {
        for (dst, src) in merged.iter_mut().zip(chunk.iter()) {
            dst.merge(src);
        }
    }

    let eff_exposure = if loaded.exposure_sec > 0.0 {
        loaded.exposure_sec
    } else if loaded.times.len() >= 2 {
        let t0 = loaded.times.iter().copied().fold(f64::INFINITY, f64::min);
        let t1 = loaded
            .times
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        (t1 - t0).max(1.0)
    } else {
        1.0
    };
    let bin_exp = eff_exposure / nbins_f;

    merged
        .into_iter()
        .zip(loaded.bands.iter())
        .map(|(acc, band)| {
            let counts = acc.bins;
            let errors: Vec<f64> = acc.bins_w2.iter().map(|&w2| w2.max(0.0).sqrt()).collect();
            let rate: Vec<f64> = counts.iter().map(|&c| c / bin_exp).collect();
            let rate_errors: Vec<f64> = errors.iter().map(|&e| e / bin_exp).collect();

            let mean_c = acc.total_weight / nbins_f;
            let (normalized, normalized_errors) = if mean_c.abs() > 1e-15 {
                (
                    counts.iter().map(|&c| c / mean_c).collect::<Vec<_>>(),
                    errors
                        .iter()
                        .map(|&e| e / mean_c.abs())
                        .collect::<Vec<_>>(),
                )
            } else {
                (vec![0.0; nbins], vec![0.0; nbins])
            };

            let (min_c, max_c, peak_bin) = counts.iter().enumerate().fold(
                (f64::INFINITY, f64::NEG_INFINITY, 0usize),
                |(mn, mx, pb), (idx, &v)| {
                    let n_mn = mn.min(v);
                    if v > mx {
                        (n_mn, v, idx)
                    } else {
                        (n_mn, mx, pb)
                    }
                },
            );
            let pulsed_fraction_minmax = if acc.num_events > 0 && (max_c + min_c).abs() > 1e-15 {
                (max_c - min_c) / (max_c + min_c)
            } else {
                0.0
            };
            let pulsed_fraction_rms = if acc.num_events > 0 {
                (normalized
                    .iter()
                    .map(|&r| (r - 1.0) * (r - 1.0))
                    .sum::<f64>()
                    / nbins_f)
                    .sqrt()
            } else {
                0.0
            };

            let z2_n = if acc.total_w2 > 0.0 {
                (2.0 / acc.total_w2)
                    * acc
                        .cos_sum
                        .iter()
                        .zip(acc.sin_sum.iter())
                        .map(|(&c, &s)| c * c + s * s)
                        .sum::<f64>()
            } else {
                0.0
            };

            let chi2 = if acc.num_events > 0 {
                if loaded.weighted {
                    counts
                        .iter()
                        .zip(acc.bins_w2.iter())
                        .map(|(&c, &var)| {
                            let d = c - mean_c;
                            if var > 0.0 {
                                d * d / var
                            } else {
                                0.0
                            }
                        })
                        .sum::<f64>()
                } else if mean_c > 0.0 {
                    counts
                        .iter()
                        .map(|&c| {
                            let d = c - mean_c;
                            d * d / mean_c
                        })
                        .sum::<f64>()
                } else {
                    0.0
                }
            } else {
                0.0
            };

            let full_label = format!("{} · {}", loaded.input_label, band.label);
            BandProfileResult {
                input_label: loaded.input_label.clone(),
                band_label: band.label.clone(),
                full_label,
                path: loaded.path.to_string_lossy().into_owned(),
                hdu_used: loaded.hdu_used.clone(),
                weighted: loaded.weighted,
                num_events: acc.num_events,
                total_weight: acc.total_weight,
                exposure_sec: eff_exposure,
                counts,
                errors,
                rate,
                rate_errors,
                normalized,
                normalized_errors,
                pulsed_fraction_minmax,
                pulsed_fraction_rms,
                z2_n,
                chi2,
                reduced_chi2: chi2 / ((nbins - 1) as f64),
                peak_phase: ((peak_bin as f64) + 0.5) / nbins_f,
                phase_lag_cycles: 0.0,
                phase_lag_ms: 0.0,
                xcorr_peak: 0.0,
            }
        })
        .collect()
}

#[inline]
fn event_matches_band(
    band: &BandSpec,
    energy_kev: Option<f64>,
    chan: Option<i32>,
    unit: Option<u8>,
    layer: Option<u8>,
) -> bool {
    if let Some(lo) = band.emin_kev {
        let Some(e) = energy_kev else { return false };
        if !e.is_finite() || e < lo {
            return false;
        }
    }
    if let Some(hi) = band.emax_kev {
        let Some(e) = energy_kev else { return false };
        if !e.is_finite() || e > hi {
            return false;
        }
    }
    if let Some(lo) = band.pi_min.or(band.chan_min) {
        let Some(c) = chan else { return false };
        if c < lo {
            return false;
        }
    }
    if let Some(hi) = band.pi_max.or(band.chan_max) {
        let Some(c) = chan else { return false };
        if c > hi {
            return false;
        }
    }
    if let Some(u_req) = band.unit {
        let Some(u) = unit else { return false };
        if u != u_req {
            return false;
        }
    }
    if let Some(l_req) = band.layer {
        let Some(l) = layer else { return false };
        if l != l_req {
            return false;
        }
    }
    true
}

/// Exact circular cross-correlation phase lag between reference profile `ref_p` and target `tgt_p`.
///
/// Positive return value `+delta_phi` (in `(-0.5, 0.5]` cycles) means `tgt_p` arrives LATER
/// in phase than `ref_p` (i.e., `tgt_p(phi) ≈ ref_p(phi - delta_phi)`).
/// Returns `(lag_cycles, normalized_peak_correlation)`.
pub fn circular_xcorr_phase_lag(ref_p: &[f64], tgt_p: &[f64]) -> (f64, f64) {
    let m = ref_p.len().min(tgt_p.len());
    if m < 2 {
        return (0.0, 0.0);
    }
    let mf = m as f64;
    let mean_a = ref_p[..m].iter().sum::<f64>() / mf;
    let mean_b = tgt_p[..m].iter().sum::<f64>() / mf;

    let a: Vec<f64> = ref_p[..m].iter().map(|&v| v - mean_a).collect();
    let b: Vec<f64> = tgt_p[..m].iter().map(|&v| v - mean_b).collect();

    let norm_a = a.iter().map(|&v| v * v).sum::<f64>();
    let norm_b = b.iter().map(|&v| v * v).sum::<f64>();
    let denom = (norm_a * norm_b).sqrt();
    if denom <= 1e-18 {
        return (0.0, 0.0);
    }

    let mut corr = vec![0.0f64; m];
    let mut best_lag = 0usize;
    let mut best_val = f64::NEG_INFINITY;

    for lag in 0..m {
        let mut s = 0.0;
        for k in 0..m {
            s += a[k] * b[(k + lag) % m];
        }
        let c = s / denom;
        corr[lag] = c;
        if c > best_val {
            best_val = c;
            best_lag = lag;
        }
    }

    // 3-point parabolic sub-bin interpolation on the circular ring
    let prev = corr[(best_lag + m - 1) % m];
    let curr = corr[best_lag];
    let next = corr[(best_lag + 1) % m];
    let curv = prev - 2.0 * curr + next;
    let delta = if curv < -1e-15 {
        (0.5 * (prev - next) / curv).clamp(-0.5, 0.5)
    } else {
        0.0
    };

    let mut lag_bins = (best_lag as f64) + delta;
    if lag_bins > 0.5 * mf {
        lag_bins -= mf;
    } else if lag_bins <= -0.5 * mf {
        lag_bins += mf;
    }
    (lag_bins / mf, best_val.clamp(-1.0, 1.0))
}

// ---------------------------------------------------------------------------
// Shared exact timing primitives
// ---------------------------------------------------------------------------

#[inline]
pub fn pulse_phase(t: f64, epoch: f64, f: f64, fdot: f64, fddot: f64) -> f64 {
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
) -> (f64, f64) {
    let mut cos_sum = vec![0.0f64; harmonics];
    let mut sin_sum = vec![0.0f64; harmonics];
    let mut bins = vec![0u64; nbins];
    let nbins_f = nbins as f64;

    for &t in times {
        let ph = pulse_phase(t, epoch, freq, fdot, fddot);
        let b = ((ph * nbins_f) as usize).min(nbins - 1);
        bins[b] += 1;

        let (s1, c1) = (TAU * ph).sin_cos();
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

    (z2_n, chi2)
}

fn fold_counts_parallel(
    times: &[f64],
    epoch: f64,
    freq: f64,
    fdot: f64,
    fddot: f64,
    nbins: usize,
) -> Vec<u64> {
    let nbins_f = nbins as f64;
    let chunks: Vec<Vec<u64>> = times
        .par_chunks(FOLD_CHUNK_SIZE)
        .map(|chunk| {
            let mut local = vec![0u64; nbins];
            for &t in chunk {
                let ph = pulse_phase(t, epoch, freq, fdot, fddot);
                let b = ((ph * nbins_f) as usize).min(nbins - 1);
                local[b] += 1;
            }
            local
        })
        .collect();
    let mut out = vec![0u64; nbins];
    for c in chunks {
        for (a, b) in out.iter_mut().zip(c) {
            *a += b;
        }
    }
    out
}

fn select_event_hdus(
    fits: &FitsFile,
    requested: Option<&str>,
    path: &Path,
) -> ApiResult<Vec<(usize, String)>> {
    if let Some(name) = requested.map(str::trim).filter(|s| !s.is_empty()) {
        if let Ok(idx) = name.parse::<usize>() {
            let h = fits
                .hdu(idx)
                .map_err(|_| ApiError::bad_request(format!("HDU {idx} not found")))?;
            let lbl = h
                .extname()
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| format!("HDU_{idx}"));
            return Ok(vec![(idx, lbl)]);
        }
        let h = fits.hdu_by_name(name).map_err(|_| {
            ApiError::bad_request(format!("HDU '{name}' not found in {}", path.display()))
        })?;
        return Ok(vec![(h.index(), name.to_string())]);
    }
    if let Ok(h) = fits.hdu_by_name("EVENTS") {
        return Ok(vec![(h.index(), "EVENTS".to_string())]);
    }
    if let Ok(h) = fits.hdu_by_name("LAXPC_EVENTS") {
        return Ok(vec![(h.index(), "LAXPC_EVENTS".to_string())]);
    }
    if let Ok(h) = fits.hdu_by_name("LAXPC_MODE_EA") {
        return Ok(vec![(h.index(), "LAXPC_MODE_EA".to_string())]);
    }
    if fits.hdu_by_name("Q0").is_ok() {
        let qs: Vec<(usize, String)> = ["Q0", "Q1", "Q2", "Q3"]
            .into_iter()
            .filter_map(|q| fits.hdu_by_name(q).ok().map(|h| (h.index(), q.to_string())))
            .collect();
        if !qs.is_empty() {
            return Ok(qs);
        }
    }
    if let Some(first_bt) = fits.iter_hdus().find(|h| {
        *h.kind() == HduKind::BinTable
            && h.as_bintable()
                .map(|bt| bt.column_index("TIME").is_some())
                .unwrap_or(false)
    }) {
        let lbl = first_bt
            .extname()
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| format!("HDU_{}", first_bt.index()));
        return Ok(vec![(first_bt.index(), lbl)]);
    }
    Err(ApiError::bad_request(
        "no BINTABLE HDU with a TIME column found in FITS file",
    ))
}

fn load_filtered_event_times(
    path: &Path,
    hdu_opt: Option<&str>,
    pi_min: Option<i32>,
    pi_max: Option<i32>,
    unit_filter: Option<u8>,
    layer_filter: Option<u8>,
) -> ApiResult<(String, Vec<f64>)> {
    let fits = open_fits_clean(path)?;
    let target_hdus = select_event_hdus(&fits, hdu_opt, path)?;
    let mut all_times = Vec::new();
    let mut labels = Vec::with_capacity(target_hdus.len());

    for (hdu_idx, label) in &target_hdus {
        let hdu = fits
            .hdu(*hdu_idx)
            .map_err(|e| ApiError::bad_request(format!("failed reading HDU '{label}': {e}")))?;
        let bt = hdu
            .as_bintable()
            .map_err(|e| ApiError::bad_request(format!("HDU '{label}' is not a BINTABLE: {e}")))?;
        let t = bt
            .read_col_f64("TIME")
            .map_err(|e| ApiError::bad_request(format!("missing TIME column in '{label}': {e}")))?;

        let chan = if pi_min.is_some() || pi_max.is_some() {
            Some(
                bt.read_col_i32("PI")
                    .or_else(|_| bt.read_col_i32("CHANNEL"))
                    .or_else(|_| bt.read_col_i32("PHA"))
                    .map_err(|_| {
                        ApiError::bad_request(format!(
                            "pi_min/pi_max requested, but HDU '{label}' has no PI/CHANNEL/PHA column"
                        ))
                    })?,
            )
        } else {
            None
        };
        let unit = if unit_filter.is_some() {
            Some(
                bt.read_col_u8("LAXPC_No")
                    .or_else(|_| bt.read_col_u8("DETECTOR_ID"))
                    .map_err(|_| {
                        ApiError::bad_request(format!(
                            "unit filter requested, but HDU '{label}' has no LAXPC_No/DETECTOR_ID column"
                        ))
                    })?,
            )
        } else {
            None
        };
        let layer = if layer_filter.is_some() {
            Some(bt.read_col_u8("Layer").map_err(|_| {
                ApiError::bad_request(format!(
                    "layer filter requested, but HDU '{label}' has no Layer column"
                ))
            })?)
        } else {
            None
        };

        for idx in 0..t.len() {
            let tv = t[idx];
            if !tv.is_finite() {
                continue;
            }
            if let (Some(cv), Some(lo)) = (&chan, pi_min) {
                if cv[idx] < lo {
                    continue;
                }
            }
            if let (Some(cv), Some(hi)) = (&chan, pi_max) {
                if cv[idx] > hi {
                    continue;
                }
            }
            if let (Some(uv), Some(u_req)) = (&unit, unit_filter) {
                if uv[idx] != u_req {
                    continue;
                }
            }
            if let (Some(lv), Some(l_req)) = (&layer, layer_filter) {
                if lv[idx] != l_req {
                    continue;
                }
            }
            all_times.push(tv);
        }
        labels.push(label.clone());
    }

    Ok((labels.join("+"), all_times))
}

fn load_gti_set(path: &Path) -> ApiResult<GtiSet> {
    let fits = open_fits_clean(path)?;
    if let Ok(hdu) = fits.hdu_by_name("GTI") {
        return GtiSet::from_hdu(&hdu).map_err(|e| {
            ApiError::bad_request(format!("failed to load GTI from {}: {e}", path.display()))
        });
    }
    for hdu in fits.iter_hdus() {
        if let Ok(bt) = hdu.as_bintable() {
            if bt.column_index("START").is_some() && bt.column_index("STOP").is_some() {
                return GtiSet::from_hdu(&hdu).map_err(|e| {
                    ApiError::bad_request(format!(
                        "failed to load GTI from {}: {e}",
                        path.display()
                    ))
                });
            }
        }
    }
    Err(ApiError::bad_request(format!(
        "no GTI HDU with START/STOP columns found in {}",
        path.display()
    )))
}

#[allow(clippy::too_many_arguments)]
fn build_frequency_grid(
    freq_hz: Option<f64>,
    period_s: Option<f64>,
    freq_min: Option<f64>,
    freq_max: Option<f64>,
    search_span_hz: f64,
    df_hz: Option<f64>,
    n_trials: Option<usize>,
    t_span_sec: f64,
) -> ApiResult<Vec<f64>> {
    let nyquist_step = 1.0 / (5.0 * t_span_sec);
    if let (Some(fmin), Some(fmax)) = (freq_min, freq_max) {
        if !fmin.is_finite() || !fmax.is_finite() || fmax <= fmin || fmin <= 0.0 {
            return Err(ApiError::bad_request(format!(
                "freq_max ({fmax}) must be finite and > freq_min ({fmin}) > 0"
            )));
        }
        let n = if let Some(nt) = n_trials {
            if nt == 0 || nt > MAX_SEARCH_TRIALS {
                return Err(ApiError::bad_request(format!(
                    "n_trials ({nt}) must be in 1..={MAX_SEARCH_TRIALS}"
                )));
            }
            nt
        } else {
            let df = df_hz.unwrap_or(nyquist_step).max(1e-12);
            let raw = ((fmax - fmin) / df).ceil();
            if !raw.is_finite() || raw >= MAX_SEARCH_TRIALS as f64 {
                MAX_SEARCH_TRIALS
            } else {
                (raw as usize).clamp(1, MAX_SEARCH_TRIALS)
            }
        };
        if n == 1 {
            return Ok(vec![0.5 * (fmin + fmax)]);
        }
        let step = (fmax - fmin) / ((n - 1) as f64);
        return Ok((0..n).map(|i| fmin + (i as f64) * step).collect());
    }

    let center_f = match (freq_hz, period_s) {
        (Some(f), _) => f,
        (None, Some(p)) => {
            if !p.is_finite() || p <= 0.0 {
                return Err(ApiError::bad_request(format!(
                    "period_s must be finite and > 0, got {p}"
                )));
            }
            1.0 / p
        }
        (None, None) => {
            return Err(ApiError::bad_request(
                "provide freq_hz, period_s, or freq_min/freq_max",
            ));
        }
    };
    if !center_f.is_finite() || center_f <= 0.0 {
        return Err(ApiError::bad_request(format!(
            "freq_hz must be finite and > 0, got {center_f}"
        )));
    }
    if let Some(nt) = n_trials {
        if nt == 0 || nt > MAX_SEARCH_TRIALS {
            return Err(ApiError::bad_request(format!(
                "n_trials ({nt}) must be in 1..={MAX_SEARCH_TRIALS}"
            )));
        }
    }

    if search_span_hz == 0.0 && n_trials.unwrap_or(1) <= 1 {
        return Ok(vec![center_f]);
    }
    let effective_span = if search_span_hz > 0.0 {
        search_span_hz
    } else {
        5.0 * nyquist_step
    };
    let fmin = (center_f - effective_span).max(1e-9);
    let fmax = center_f + effective_span;
    if !fmin.is_finite() || !fmax.is_finite() || fmax <= fmin {
        return Err(ApiError::bad_request(format!(
            "invalid frequency range [{fmin}, {fmax}]"
        )));
    }
    let n = if let Some(nt) = n_trials {
        nt
    } else {
        let df = df_hz.unwrap_or(nyquist_step).max(1e-12);
        let raw = ((fmax - fmin) / df).ceil();
        if !raw.is_finite() || raw >= MAX_SEARCH_TRIALS as f64 {
            MAX_SEARCH_TRIALS
        } else {
            (raw as usize).clamp(3, MAX_SEARCH_TRIALS)
        }
    };
    if n == 1 {
        return Ok(vec![center_f]);
    }
    let step = (fmax - fmin) / ((n - 1) as f64);
    Ok((0..n).map(|i| fmin + (i as f64) * step).collect())
}

/// Exact log-space survival probability P(Z^2_n >= z) = exp(-z/2) * sum_{j=0}^{n-1} (z/2)^j / j!.
pub fn z2_fap(z2: f64, harmonics: usize) -> f64 {
    if !z2.is_finite() || z2 <= 0.0 || harmonics == 0 {
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

/// Convert tail probability `p` in `(0, 1)` to one-sided Gaussian sigma.
pub fn fap_to_sigma(p: f64) -> f64 {
    if !p.is_finite() || p >= 0.5 {
        return 0.0;
    }
    let p = p.max(1e-300);
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
