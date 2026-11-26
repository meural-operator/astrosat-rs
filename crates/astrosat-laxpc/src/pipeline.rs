//! Stage 7: End-to-end `rayon`-parallelized AstroSat LAXPC Level-1 -> Level-2 pipeline orchestrator.
//!
//! Discovers single- or multi-orbit AstroSat Level-1 LAXPC directories (or accepts explicit orbit
//! file sets), processes all requested LAXPC units (LAXPC10/20/30) in parallel, and writes
//! OGIP-compliant Level-2 Format B and optional Format A products + `AstroSat.orb`.

use crate::background::{
    apply_lightcurve_gain_shift, background_lightcurve_shift_fraction, compute_gain_shift,
    shift_background_spectrum, BackgroundCatalog, BackgroundModel2D, BackgroundPha, ChannelMap,
    GainShiftSolution, GainTable, HvJumpTable,
};
use crate::bb::BbcSeries;
use crate::ea::decode_ea_orbits;
use crate::error::{LaxpcError, Result};
use crate::frames::{sort_and_dedup_frames, FrameDiagnostics, FrameMode, OrderedFrame};
use crate::gti::{
    generate_standard_gti, interp_lon_deg, load_laxpc_tct_files, load_merged_mkf_files, utc_to_ins,
    GtiConfig, GtiResult,
};
use crate::products::{
    bin_events_raw, build_calibrated_events, build_lightcurve, build_spectrum, write_ascii_lc,
    write_ascii_pha, write_format_a_event_fits_with_channels, write_format_b_event_fits,
    write_gti_fits, write_lightcurve_fits, write_orbit_fits, write_spectrum_fits, CalibratedEvent,
    LaxpcLightCurve, LaxpcSpectrum,
};
use crate::types::{
    utc_to_mission_day, EventFilter, LaxpcUnit, DEFAULT_DEADTIME_SEC, L1_FRAME_SIZE,
};
use astrosat_core::gti::{GtiInterval, GtiSet};
use astrosat_core::mkf::MkfTable;
use astrosat_fits::{FitsFile, HduKind, OgipMeta};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

/// File paths for a single orbit of one LAXPC unit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitInputFiles {
    /// Orbit identifier (e.g. `"01509"`).
    pub orbit_id: String,
    /// Level-1 `.mkf` housekeeping file.
    pub mkf_path: PathBuf,
    /// Level-1 `.tct` Time Correlation Table file.
    pub tct_path: PathBuf,
    /// Level-1 `modeEA` FITS file (`*lxp{1,2,3}EA_level1.fits`).
    pub ea_path: PathBuf,
    /// Optional Level-1 `modeBB` FITS file (`*lxp{1,2,3}BB_level1.fits`).
    pub bb_path: Option<PathBuf>,
}

/// Complete configuration for running the LAXPC Level-1 -> Level-2 pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    /// Optional Level-1 root directory for automatic orbit/unit discovery.
    pub input_dir: Option<PathBuf>,
    /// Optional explicit per-unit orbit file lists (overrides `input_dir` discovery for listed units).
    pub orbits_by_unit: BTreeMap<LaxpcUnit, Vec<OrbitInputFiles>>,
    /// Output directory where Level-2 FITS/ASCII products are written.
    pub output_dir: PathBuf,
    /// Which LAXPC units to process (default: `[Lxp10, Lxp20, Lxp30]`).
    pub units: Vec<LaxpcUnit>,
    /// Light curve time bin width in seconds (default `1.0`).
    pub tbin: f64,
    /// Non-paralyzable deadtime constant in seconds (default `43e-6`, configurable `42e-6` / `43e-6`).
    pub deadtime_sec: f64,
    /// Spectrum channel count (`256`, `512`, `1024`); `None` uses unit default (`512` for Lxp10/30, `256` for Lxp20).
    pub spectrum_channels: Option<usize>,
    /// Event/anode/energy selection filter.
    pub filter: EventFilter,
    /// Standard GTI screening parameters.
    pub gti_config: GtiConfig,
    /// Optional user GTI in UTC seconds; if provided, intersected with standard GTI (or used directly if `use_user_gti_only`).
    pub user_gti_utc: Option<GtiSet>,
    /// If `true` and `user_gti_utc` is `Some`, bypass auto GTI and use `user_gti_utc` directly.
    pub use_user_gti_only: bool,
    /// Optional path to LAXPC CALDB / reference directory (`gaina8.dat`, `lxhv.dat`, `lxp{10,20,30}.channel`, `backfit.*`).
    pub caldb_dir: Option<PathBuf>,
    /// Optional explicit 2D background model (`backfit.lxp?0*`) file.
    pub background_model_path: Option<PathBuf>,
    /// Optional explicit background `.pha` file (`backlxp?0*.pha`) for background spectrum generation & shifting.
    pub background_pha_path: Option<PathBuf>,
    /// Optional explicit 30 keV calibration peak channel (`chb`) for a custom background `.pha`
    /// (overrides automatic `backshiftv3.inp` catalog lookup).
    pub background_peak_ch: Option<f64>,
    /// Write per-unit Format B `.evn` event FITS files (default `true`).
    pub write_format_b_events: bool,
    /// Write combined multi-unit Format A `level2.event.fits` (`TIME`, `Channel`, `Layer`, `LAXPC_No`) (default `false`).
    pub write_format_a_combined: bool,
    /// Write `AstroSat.orb` orbit FITS file interpolated from MKF (default `true`).
    pub write_orbit_file: bool,
    /// Also emit ASCII `.pha` and `.lc` diagnostics alongside FITS products (default `false`).
    pub write_ascii_products: bool,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            input_dir: None,
            orbits_by_unit: BTreeMap::new(),
            output_dir: PathBuf::from("laxpc_l2_out"),
            units: LaxpcUnit::ALL.to_vec(),
            tbin: 1.0,
            deadtime_sec: DEFAULT_DEADTIME_SEC,
            spectrum_channels: None,
            filter: EventFilter::default(),
            gti_config: GtiConfig::default(),
            user_gti_utc: None,
            use_user_gti_only: false,
            caldb_dir: None,
            background_model_path: None,
            background_pha_path: None,
            background_peak_ch: None,
            write_format_b_events: true,
            write_format_a_combined: false,
            write_orbit_file: true,
            write_ascii_products: false,
        }
    }
}

/// Summary and output file paths produced for one LAXPC unit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnitPipelineResult {
    pub unit: LaxpcUnit,
    pub num_orbits: usize,
    pub ea_frame_diagnostics: FrameDiagnostics,
    pub bb_frame_diagnostics: Option<FrameDiagnostics>,
    pub total_valid_events: usize,
    pub gti_filtered_events: usize,
    pub gti: GtiResult,
    pub spectrum: LaxpcSpectrum,
    pub lightcurve_bins: usize,
    pub lightcurve_mean_rate_cps: f64,
    pub gain_shift: Option<GainShiftSolution>,
    pub recommended_rmf: Option<String>,
    pub event_fits_path: Option<PathBuf>,
    pub spectrum_fits_path: PathBuf,
    pub lightcurve_fits_path: PathBuf,
    pub gti_fits_path: PathBuf,
    pub bkg_spectrum_fits_path: Option<PathBuf>,
    pub bkg_lightcurve_fits_path: Option<PathBuf>,
}

/// Aggregate result of running the multi-unit LAXPC pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiUnitPipelineResult {
    pub unit_results: Vec<UnitPipelineResult>,
    pub orbit_fits_path: Option<PathBuf>,
    pub format_a_event_fits_path: Option<PathBuf>,
    pub format_a_total_events: usize,
}

impl MultiUnitPipelineResult {
    /// Serialize the pipeline summary to pretty-printed JSON (useful for CLI `--json` output).
    pub fn to_json_pretty(&self) -> Result<String> {
        serde_json::to_string_pretty(self).map_err(|e| LaxpcError::Config(e.to_string()))
    }
}

/// Automatically discover AstroSat LAXPC Level-1 orbit files for a specific `unit` under `root`.
///
/// Works whether `root` points to the archive root (`20160107_..._level1`), `.../laxpc`,
/// or a single orbit directory (`.../01509`).
pub fn discover_l1_orbits(root: &Path, unit: LaxpcUnit) -> Result<Vec<OrbitInputFiles>> {
    if !root.exists() {
        return Err(LaxpcError::MissingInput(root.to_path_buf()));
    }

    let mut orbit_dirs = Vec::new();
    collect_orbit_dirs(root, &mut orbit_dirs)?;
    orbit_dirs.sort();
    orbit_dirs.dedup();

    let prefix = unit.prefix();
    let mut orbits = Vec::new();

    for dir in orbit_dirs {
        let mkf = find_file_matching(&dir, |p| {
            p.extension().and_then(|e| e.to_str()) == Some("mkf")
        });
        let aux_dir = dir.join("aux");
        let tct = if aux_dir.is_dir() {
            find_file_matching(&aux_dir, |p| {
                p.extension().and_then(|e| e.to_str()) == Some("tct")
            })
        } else {
            find_file_matching(&dir, |p| {
                p.extension().and_then(|e| e.to_str()) == Some("tct")
            })
        };

        let ea_dir = dir.join(prefix).join("modeEA");
        let ea = if ea_dir.is_dir() {
            find_file_matching(&ea_dir, |p| {
                let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
                name.ends_with(".fits") && !name.starts_with('.')
            })
        } else {
            None
        };

        let bb_dir = dir.join(prefix).join("modeBB");
        let bb = if bb_dir.is_dir() {
            find_file_matching(&bb_dir, |p| {
                let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
                name.ends_with(".fits") && !name.starts_with('.')
            })
        } else {
            None
        };

        if let (Some(mkf_path), Some(tct_path), Some(ea_path)) = (mkf, tct, ea) {
            let orbit_id = dir
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("orbit")
                .to_string();
            orbits.push(OrbitInputFiles {
                orbit_id,
                mkf_path,
                tct_path,
                ea_path,
                bb_path: bb,
            });
        }
    }

    if orbits.is_empty() {
        return Err(LaxpcError::MissingInput(root.join(prefix)));
    }
    Ok(orbits)
}

fn collect_orbit_dirs(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    // An orbit directory directly contains a `.mkf` file and at least one `lxp{1,2,3}` subdirectory
    let has_lxp = dir.join("lxp1").is_dir() || dir.join("lxp2").is_dir() || dir.join("lxp3").is_dir();
    if has_lxp {
        out.push(dir.to_path_buf());
        return Ok(());
    }
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            collect_orbit_dirs(&path, out)?;
        }
    }
    Ok(())
}

fn find_file_matching<F>(dir: &Path, pred: F) -> Option<PathBuf>
where
    F: Fn(&Path) -> bool,
{
    let mut matches: Vec<PathBuf> = fs::read_dir(dir)
        .ok()?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.is_file() && pred(p))
        .collect();
    matches.sort();
    matches.into_iter().next()
}

/// Read all 2242-byte raw frames from a Level-1 EA or BB FITS BINTABLE HDU (`NAXIS1 = 2242`)
/// and sort/deduplicate them.
pub fn read_and_sort_l1_fits_frames(
    path: &Path,
    unit: LaxpcUnit,
    mode: FrameMode,
) -> Result<(Vec<OrderedFrame>, FrameDiagnostics)> {
    let fits = FitsFile::open(path).map_err(|e| LaxpcError::Fits(e.to_string()))?;
    let ext_preferred = match mode {
        FrameMode::Ea => "LAXPC_MODEEA",
        FrameMode::Bb => "LAXPC_MODEBB",
    };
    let hdu = fits
        .hdu_by_name(ext_preferred)
        .or_else(|_| {
            fits.iter_hdus()
                .find(|h| *h.kind() == HduKind::BinTable)
                .ok_or_else(|| astrosat_fits::FitsError::HduNotFound(ext_preferred.into()))
        })
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    let bt = hdu
        .as_bintable()
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;

    if bt.row_len() < L1_FRAME_SIZE {
        return Err(LaxpcError::Frame(format!(
            "{}: BINTABLE NAXIS1={} < {L1_FRAME_SIZE}",
            path.display(),
            bt.row_len()
        )));
    }

    let nrows = bt.nrows();
    let raw_iter = (0..nrows).map(|r| bt.row_bytes(r));
    sort_and_dedup_frames(raw_iter, Some(unit), Some(mode))
}

/// Extract header metadata (`OBJECT`, `OBS_ID`, `RA_NOM`, `DEC_NOM`, `DATE-OBS`, `DATE-END`)
/// from the first/last MKF files.
fn read_mkf_header_meta(mkf_paths: &[PathBuf]) -> OgipMeta {
    let mut meta = OgipMeta {
        instrume: "LAXPC".to_string(),
        creator: "astrosat-rs/astrosat-laxpc 0.1.0".to_string(),
        origin: "ASTROSAT-LAXPC-L2".to_string(),
        ..Default::default()
    };
    if let Some(first) = mkf_paths.first() {
        if let Ok(fits) = FitsFile::open(first) {
            if let Ok(h0) = fits.hdu(0) {
                let hdr = h0.header();
                meta.object = hdr.get_str("OBJECT").map(|s| s.trim().to_string());
                meta.obs_id = hdr.get_str("OBS_ID").map(|s| s.trim().to_string());
                meta.ra_obj = hdr.get_f64("RA_NOM").or_else(|| hdr.get_f64("RA_OBJ"));
                meta.dec_obj = hdr.get_f64("DEC_NOM").or_else(|| hdr.get_f64("DEC_OBJ"));
                meta.date_obs = hdr.get_str("DATE-OBS").map(|s| s.trim().to_string());
                meta.date_end = hdr.get_str("DATE-END").map(|s| s.trim().to_string());
            }
        }
    }
    if let Some(last) = mkf_paths.last() {
        if let Ok(fits) = FitsFile::open(last) {
            if let Ok(h0) = fits.hdu(0) {
                if let Some(de) = h0.header().get_str("DATE-END") {
                    meta.date_end = Some(de.trim().to_string());
                }
            }
        }
    }
    meta
}

/// Internal per-unit execution returning both public summary and calibrated events (for Format A merge).
fn run_unit_internal(
    unit: LaxpcUnit,
    orbits: &[OrbitInputFiles],
    config: &PipelineConfig,
) -> Result<(UnitPipelineResult, Vec<CalibratedEvent>, MkfTable, OgipMeta)> {
    if orbits.is_empty() {
        return Err(LaxpcError::Config(format!("no orbits provided for {unit}")));
    }

    fs::create_dir_all(&config.output_dir)?;

    let mkf_paths: Vec<PathBuf> = orbits.iter().map(|o| o.mkf_path.clone()).collect();
    let tct_paths: Vec<PathBuf> = orbits.iter().map(|o| o.tct_path.clone()).collect();
    let mut base_meta = read_mkf_header_meta(&mkf_paths);
    base_meta.detnam = Some(unit.detnam().to_string());

    // Resolve optional channel map for energy->channel selection and gain shift
    let channel_map = resolve_channel_map(config.caldb_dir.as_deref(), unit);
    let mut effective_filter = config.filter.clone();
    if let Some((e_lo, e_hi)) = config.filter.energy_range_kev {
        let cmap = channel_map.as_ref().ok_or_else(|| {
            LaxpcError::Config(format!(
                "energy_range_kev ({e_lo}..={e_hi} keV) requested for {unit}, but {}0.channel was not found in caldb_dir ({:?})",
                unit.prefix(),
                config.caldb_dir
            ))
        })?;
        let (p_lo, p_hi) = cmap.energy_to_pha_1024(e_lo, e_hi);
        effective_filter.pha_min = p_lo;
        effective_filter.pha_max = p_hi;
    }

    // Parallel read & frame-sort across orbits for EA and BB
    let ea_batches_res: Result<Vec<(Vec<OrderedFrame>, FrameDiagnostics)>> = orbits
        .par_iter()
        .map(|o| read_and_sort_l1_fits_frames(&o.ea_path, unit, FrameMode::Ea))
        .collect();
    let ea_batches = ea_batches_res?;

    let mut ea_diag = FrameDiagnostics::default();
    let mut ea_orbit_frames: Vec<Vec<OrderedFrame>> = Vec::with_capacity(ea_batches.len());
    for (frames, d) in ea_batches {
        accumulate_diag(&mut ea_diag, d);
        ea_orbit_frames.push(frames);
    }

    let bb_paths: Vec<&Path> = orbits
        .iter()
        .filter_map(|o| o.bb_path.as_deref())
        .collect();
    let (bbc_series, bb_diag) = if !bb_paths.is_empty() {
        let bb_batches_res: Result<Vec<(Vec<OrderedFrame>, FrameDiagnostics)>> = bb_paths
            .par_iter()
            .map(|p| read_and_sort_l1_fits_frames(p, unit, FrameMode::Bb))
            .collect();
        let bb_batches = bb_batches_res?;
        let mut diag = FrameDiagnostics::default();
        let mut bb_orbit_frames = Vec::with_capacity(bb_batches.len());
        for (frames, d) in bb_batches {
            accumulate_diag(&mut diag, d);
            bb_orbit_frames.push(frames);
        }
        (Some(BbcSeries::from_ordered_orbits(&bb_orbit_frames)), Some(diag))
    } else {
        (None, None)
    };

    // Load TCT and MKF
    let tct = load_laxpc_tct_files(&tct_paths)?;
    let mkf = load_merged_mkf_files(&mkf_paths)?;

    // Stage 2: Decode EA events
    let ea_stream = decode_ea_orbits(&ea_orbit_frames, unit);
    tracing::info!(
        unit = %unit,
        retained_frames = ea_diag.retained_frames,
        missing_frames = ea_diag.total_missing_frames,
        valid_events = ea_stream.total_valid_events,
        "decoded LAXPC EA frames"
    );

    // Stage 4: Raw 1-s (or tbin) binning + Standard GTI generation
    let gti_tbin = if config.tbin < 0.2 { 1.0 } else { config.tbin };
    let gti_raw_bins = bin_events_raw(&ea_stream.events, gti_tbin, &effective_filter)?;

    let mut gti_cfg = config.gti_config.clone();
    gti_cfg.deadtime_sec = config.deadtime_sec;

    let gti_res = if config.use_user_gti_only {
        if let Some(ugti) = &config.user_gti_utc {
            let gti_ins = GtiSet::new(
                ugti.intervals()
                    .iter()
                    .map(|iv| GtiInterval::new(utc_to_ins(&tct, iv.start), utc_to_ins(&tct, iv.stop)))
                    .collect(),
            );
            let pairs_ins: Vec<(f64, f64)> = gti_ins
                .intervals()
                .iter()
                .map(|iv| (iv.start, iv.stop))
                .collect();
            let gross = gti_ins.exposure();
            let gap_loss = crate::products::compute_gti_frame_gap_loss(&pairs_ins, &ea_stream.frame_gaps);
            GtiResult {
                gti_ins,
                gti_utc: ugti.clone(),
                gross_exposure_sec: gross,
                frame_gap_loss_sec: gap_loss,
                net_exposure_sec: (gross - gap_loss).max(0.0),
                uld_quartiles: Default::default(),
                max_ang_offset_deg: 0.0,
                max_ang_offset_time_utc: 0.0,
                tct_rejected_glitches: tct.rejected_glitch_count(),
            }
        } else {
            generate_standard_gti(
                unit,
                &gti_raw_bins,
                &tct,
                &mkf,
                bbc_series.as_ref(),
                &ea_stream.frame_gaps,
                &gti_cfg,
                None,
            )?
        }
    } else {
        generate_standard_gti(
            unit,
            &gti_raw_bins,
            &tct,
            &mkf,
            bbc_series.as_ref(),
            &ea_stream.frame_gaps,
            &gti_cfg,
            config.user_gti_utc.as_ref(),
        )?
    };

    let gti_ins_pairs: Vec<(f64, f64)> = gti_res
        .gti_ins
        .intervals()
        .iter()
        .map(|iv| (iv.start, iv.stop))
        .collect();
    let gti_utc_pairs: Vec<(f64, f64)> = gti_res
        .gti_utc
        .intervals()
        .iter()
        .map(|iv| (iv.start, iv.stop))
        .collect();

    let tstart_utc = gti_res
        .gti_utc
        .tstart()
        .unwrap_or_else(|| tct.interpolate(ea_stream.t_ins_min));
    let tstop_utc = gti_res
        .gti_utc
        .tstop()
        .unwrap_or_else(|| tct.interpolate(ea_stream.t_ins_max));
    let tmid_utc = 0.5 * (tstart_utc + tstop_utc);
    let mission_day = utc_to_mission_day(tmid_utc).max(0.0);

    base_meta.tstart = tstart_utc;
    base_meta.tstop = tstop_utc;
    base_meta.exposure = gti_res.net_exposure_sec;
    base_meta.ontime = Some(gti_res.gross_exposure_sec);
    base_meta.livetime = Some(gti_res.net_exposure_sec);

    // Stage 5: Spectrum & Light Curve generation
    let num_channels = config
        .spectrum_channels
        .unwrap_or_else(|| unit.default_spectrum_channels());

    let spectrum = build_spectrum(
        &ea_stream.events,
        &gti_ins_pairs,
        &ea_stream.frame_gaps,
        unit,
        &effective_filter,
        num_channels,
        config.deadtime_sec,
    )?;

    let lc_raw_bins = if (config.tbin - gti_tbin).abs() < 1e-12 {
        gti_raw_bins
    } else {
        bin_events_raw(&ea_stream.events, config.tbin, &effective_filter)?
    };

    // Optional Stage 6: 2D background model evaluation per bin + gain alignment + RMF recommendation
    let bkg_model = resolve_background_model(config, unit);
    let raw_bkg_per_bin: Option<Vec<f64>> = bkg_model.as_ref().map(|model| {
        let mkf_t = mkf.times();
        let lon_col = mkf.column_any(&["EARTHLON", "SAT_LON", "LONGITUDE", "LON"]);
        let lat_col = mkf.column_any(&["EARTHLAT", "SAT_LAT", "LATITUDE", "LAT"]);
        (0..lc_raw_bins.all_main_counts.len())
            .into_par_iter()
            .map(|i| {
                let tc = lc_raw_bins.t0_ins + (i as f64 + 0.5) * lc_raw_bins.tbin;
                let xb = tct.interpolate(tc);
                let mut lon = lon_col
                    .map(|c| interp_lon_deg(mkf_t, c, xb))
                    .unwrap_or(100.0);
                if lon < 0.0 {
                    lon += 360.0;
                }
                let lat = lat_col
                    .map(|c| crate::bb::interp_linear(mkf_t, c, xb))
                    .unwrap_or(0.0);
                model.evaluate(lon, lat, xb)
            })
            .collect()
    });

    let mut lightcurve = build_lightcurve(
        &lc_raw_bins,
        &gti_ins_pairs,
        &ea_stream.frame_gaps,
        unit,
        config.deadtime_sec,
        |t_ins| tct.interpolate(t_ins),
        raw_bkg_per_bin.as_deref(),
    );

    // Gain table + HV jump table -> RMF recommendation & optional background shifting
    let (gain_shift, recommended_rmf, bkg_spec_path, bkg_lc_path) =
        run_optional_background_and_gain(
            unit,
            mission_day,
            config,
            &effective_filter,
            num_channels,
            channel_map.as_ref(),
            &mut lightcurve,
            &gti_utc_pairs,
            &base_meta,
        )?;

    // Calibrated UTC events for Format B (unit-native channels) and Format A (uniform channel scale)
    let format_b_events = if config.write_format_b_events {
        build_calibrated_events(
            &ea_stream.events,
            &gti_ins_pairs,
            &effective_filter,
            num_channels,
            |t_ins| tct.interpolate(t_ins),
        )
    } else {
        Vec::new()
    };

    let fa_channels = config.spectrum_channels.unwrap_or(1024);
    let format_a_events = if config.write_format_a_combined {
        if config.write_format_b_events && fa_channels == num_channels {
            format_b_events.clone()
        } else {
            build_calibrated_events(
                &ea_stream.events,
                &gti_ins_pairs,
                &effective_filter,
                fa_channels,
                |t_ins| tct.interpolate(t_ins),
            )
        }
    } else {
        Vec::new()
    };
    let gti_filtered_count = if config.write_format_b_events {
        format_b_events.len()
    } else {
        format_a_events.len()
    };

    // Write Format B per-unit FITS products
    let prefix = unit.prefix();
    let respfile = recommended_rmf.as_deref().unwrap_or("NONE");
    let backfile = bkg_spec_path
        .as_ref()
        .and_then(|p| p.file_name())
        .and_then(|n| n.to_str())
        .unwrap_or("NONE");

    let spec_path = config.output_dir.join(format!("{prefix}level2.spec"));
    write_spectrum_fits(
        &spec_path,
        unit,
        &spectrum.selected_rate,
        &spectrum.selected_stat_err,
        spectrum.exposure_sec,
        config.deadtime_sec,
        &effective_filter,
        respfile,
        backfile,
        &base_meta,
    )?;

    let lc_path = config.output_dir.join(format!("{prefix}level2.lcurv"));
    write_lightcurve_fits(
        &lc_path,
        &lightcurve,
        &lightcurve.rate_cps,
        &lightcurve.error_cps,
        &gti_utc_pairs,
        &effective_filter,
        config.deadtime_sec,
        false,
        &base_meta,
    )?;

    let gti_path = config.output_dir.join(format!("{prefix}level2.gtis"));
    write_gti_fits(&gti_path, &gti_utc_pairs, &base_meta)?;

    let event_path = if config.write_format_b_events {
        let p = config.output_dir.join(format!("{prefix}level2.evn"));
        write_format_b_event_fits(
            &p,
            &format_b_events,
            &gti_utc_pairs,
            unit,
            &effective_filter,
            num_channels,
            &base_meta,
        )?;
        Some(p)
    } else {
        None
    };

    if config.write_ascii_products {
        let pha_path = config.output_dir.join(format!("{prefix}level2.pha"));
        write_ascii_pha(&pha_path, &spectrum, mission_day.round() as i32, effective_filter.single_only)?;
        let asc_lc_path = config.output_dir.join(format!("{prefix}level2.lc"));
        write_ascii_lc(&asc_lc_path, &lightcurve)?;
    }

    let result = UnitPipelineResult {
        unit,
        num_orbits: orbits.len(),
        ea_frame_diagnostics: ea_diag,
        bb_frame_diagnostics: bb_diag,
        total_valid_events: ea_stream.total_valid_events,
        gti_filtered_events: gti_filtered_count,
        gti: gti_res,
        spectrum,
        lightcurve_bins: lightcurve.time_utc.len(),
        lightcurve_mean_rate_cps: lightcurve.mean_rate_cps,
        gain_shift,
        recommended_rmf,
        event_fits_path: event_path,
        spectrum_fits_path: spec_path,
        lightcurve_fits_path: lc_path,
        gti_fits_path: gti_path,
        bkg_spectrum_fits_path: bkg_spec_path,
        bkg_lightcurve_fits_path: bkg_lc_path,
    };

    Ok((result, format_a_events, mkf, base_meta))
}

type OptionalGainOutput = (
    Option<GainShiftSolution>,
    Option<String>,
    Option<PathBuf>,
    Option<PathBuf>,
);

type InternalUnitRun = (UnitPipelineResult, Vec<CalibratedEvent>, MkfTable, OgipMeta);

#[allow(clippy::too_many_arguments)]
fn run_optional_background_and_gain(
    unit: LaxpcUnit,
    mission_day: f64,
    config: &PipelineConfig,
    filter: &EventFilter,
    num_channels: usize,
    channel_map: Option<&ChannelMap>,
    lightcurve: &mut LaxpcLightCurve,
    gti_utc_pairs: &[(f64, f64)],
    base_meta: &OgipMeta,
) -> Result<OptionalGainOutput> {
    let gain_table = config
        .caldb_dir
        .as_deref()
        .map(|d| d.join("gaina8.dat"))
        .filter(|p| p.is_file())
        .map(|p| GainTable::from_file(&p))
        .transpose()?;

    let hv_table = config
        .caldb_dir
        .as_deref()
        .and_then(|d| {
            let p1 = d.join("lxhv.dat");
            if p1.is_file() {
                Some(p1)
            } else {
                let p2 = d.join("lxhv.data");
                p2.is_file().then_some(p2)
            }
        })
        .map(|p| HvJumpTable::from_file(&p))
        .transpose()?;

    let resolved_bkg_pha = resolve_background_pha(config, unit)?;
    let (bkg_pha_file, bkg_pha) = match resolved_bkg_pha {
        Some((p, pha)) => (Some(p), Some(pha)),
        None => (None, None),
    };

    let gain_shift = if let Some(gt) = &gain_table {
        let bkg_day = bkg_pha
            .as_ref()
            .map(|b| b.mission_day)
            .unwrap_or(mission_day);
        let bkg_peak_ch = config.background_peak_ch.or_else(|| {
            if bkg_pha.is_some() {
                let catalog = BackgroundCatalog::load_or_builtin(&[
                    config.caldb_dir.as_deref(),
                    bkg_pha_file.as_deref().and_then(Path::parent),
                ]);
                catalog.lookup_peak_30kev(unit, bkg_pha_file.as_deref(), bkg_day)
            } else {
                None
            }
        });
        compute_gain_shift(
            unit,
            mission_day,
            bkg_day,
            bkg_peak_ch,
            gt,
            hv_table.as_ref(),
            filter.anodes.is_layer1_only(),
            filter.single_only,
        )
        .ok()
    } else {
        None
    };

    let recommended_rmf = gain_shift.as_ref().map(|g| g.recommended_rmf.clone());

    let mut bkg_spec_path = None;
    let mut bkg_lc_path = None;

    let epoch_channel_map = config
        .caldb_dir
        .as_deref()
        .and_then(|d| ChannelMap::from_caldb_for_day(d, unit, mission_day));
    let effective_cmap = epoch_channel_map.as_ref().or(channel_map);

    if let (Some(pha), Some(gs), Some(cmap)) = (&bkg_pha, &gain_shift, effective_cmap) {
        if pha.rates.len() == num_channels {
            // Scale background by mean modeled rate vs background PHA total rate when 2D model is available
            let orig_total: f64 = pha.rates.iter().map(|r| r[10]).sum();
            let scale_2d = lightcurve
                .bkg_rate_cps
                .as_ref()
                .filter(|v| !v.is_empty() && orig_total > 0.0)
                .map(|v| {
                    let mean_bk = v.iter().sum::<f64>() / (v.len() as f64);
                    mean_bk / orig_total
                })
                .unwrap_or(1.0);

            let scaled_orig: Vec<[f64; 11]> = pha
                .rates
                .iter()
                .map(|row| {
                    let mut out = [0.0; 11];
                    for k in 0..11 {
                        out[k] = row[k] * scale_2d;
                    }
                    out
                })
                .collect();

            let shifted = shift_background_spectrum(
                unit,
                gs.delta_ch_30kev,
                gs.secular_scale,
                mission_day,
                cmap,
                &scaled_orig,
            )?;

            let ncf = (1024 / num_channels).max(1);
            let ch_lo = (filter.pha_min as usize) / ncf;
            let ch_hi = (filter.pha_max as usize) / ncf;
            let (_cdel, fshbk) = background_lightcurve_shift_fraction(
                &scaled_orig,
                &shifted,
                &filter.anodes,
                ch_lo,
                ch_hi,
            );

            let mut bkg_sel_rate = Vec::with_capacity(num_channels);
            let mut bkg_sel_err = Vec::with_capacity(num_channels);
            for row in &shifted {
                let s: f64 = (1..=10u8)
                    .filter(|&a| filter.anodes.matches(a))
                    .map(|a| row[(a - 1) as usize])
                    .sum();
                bkg_sel_rate.push(s);
                bkg_sel_err.push((s / pha.exposure.max(1.0)).max(0.0).sqrt());
            }

            let out_p = config
                .output_dir
                .join(format!("{}level2back_shifted.spec", unit.prefix()));
            let mut bk_meta = base_meta.clone();
            bk_meta.object = Some("BACKGROUND".to_string());
            write_spectrum_fits(
                &out_p,
                unit,
                &bkg_sel_rate,
                &bkg_sel_err,
                pha.exposure,
                config.deadtime_sec,
                filter,
                "NONE",
                "NONE",
                &bk_meta,
            )?;
            bkg_spec_path = Some(out_p);

            // Apply gain-shift adjustment to modeled background light curve and write `.lcbk`
            if let Some(mut bk_lc) = lightcurve.bkg_rate_cps.take() {
                let sel_frac = if orig_total > 0.0 {
                    let sbi: f64 = (ch_lo..=ch_hi.min(num_channels - 1))
                        .map(|i| {
                            (1..=10u8)
                                .filter(|&a| filter.anodes.matches(a))
                                .map(|a| pha.rates[i][(a - 1) as usize])
                                .sum::<f64>()
                        })
                        .sum();
                    sbi / orig_total
                } else {
                    1.0
                };
                for r in bk_lc.iter_mut() {
                    *r *= sel_frac;
                }
                apply_lightcurve_gain_shift(&mut bk_lc, fshbk);
                let bk_err: Vec<f64> = bk_lc.iter().map(|&r| (r / 32.0).max(0.0).sqrt()).collect();
                let lcbk_p = config
                    .output_dir
                    .join(format!("{}level2.lcbk", unit.prefix()));
                write_lightcurve_fits(
                    &lcbk_p,
                    lightcurve,
                    &bk_lc,
                    &bk_err,
                    gti_utc_pairs,
                    filter,
                    config.deadtime_sec,
                    false,
                    base_meta,
                )?;
                lightcurve.bkg_rate_cps = Some(bk_lc);
                bkg_lc_path = Some(lcbk_p);
            }
        }
    }

    // If a 2D background model was evaluated (`lightcurve.bkg_rate_cps` is `Some`) but `.lcbk`
    // was not written above (e.g. `background_pha_path` was omitted), still write `lxp?level2.lcbk`
    // directly from the 2D model (`laxpcl1.f:2095-2102`).
    if bkg_lc_path.is_none() {
        if let Some(bk_lc) = lightcurve.bkg_rate_cps.as_ref() {
            let bk_err: Vec<f64> = bk_lc.iter().map(|&r| (r / 32.0).max(0.0).sqrt()).collect();
            let lcbk_p = config
                .output_dir
                .join(format!("{}level2.lcbk", unit.prefix()));
            write_lightcurve_fits(
                &lcbk_p,
                lightcurve,
                bk_lc,
                &bk_err,
                gti_utc_pairs,
                filter,
                config.deadtime_sec,
                false,
                base_meta,
            )?;
            bkg_lc_path = Some(lcbk_p);
        }
    }

    Ok((gain_shift, recommended_rmf, bkg_spec_path, bkg_lc_path))
}

fn resolve_background_pha(
    config: &PipelineConfig,
    unit: LaxpcUnit,
) -> Result<Option<(PathBuf, BackgroundPha)>> {
    let Some(p) = &config.background_pha_path else {
        return Ok(None);
    };
    // Auto-substitute `lxp10` / `lxp1` -> target unit (`laxpcl1.f:1979-1981`) so multi-unit runs
    // never load LAXPC10's background spectrum for LAXPC20/LAXPC30.
    let s = p.to_string_lossy();
    let candidate = if unit != LaxpcUnit::Lxp10 {
        if s.contains("lxp10") {
            PathBuf::from(s.replace("lxp10", &format!("{}0", unit.prefix())))
        } else if s.contains("lxp1") {
            PathBuf::from(s.replace("lxp1", unit.prefix()))
        } else {
            p.clone()
        }
    } else {
        p.clone()
    };
    if candidate.is_file() {
        Ok(Some((candidate.clone(), BackgroundPha::from_file(&candidate)?)))
    } else {
        Ok(None)
    }
}

fn resolve_channel_map(caldb_dir: Option<&Path>, unit: LaxpcUnit) -> Option<ChannelMap> {
    let dir = caldb_dir?;
    let name = format!("{}0.channel", unit.prefix()); // lxp10.channel, lxp20.channel, lxp30.channel
    let p = dir.join(name);
    if p.is_file() {
        ChannelMap::from_file(&p).ok()
    } else {
        None
    }
}

fn resolve_background_model(config: &PipelineConfig, unit: LaxpcUnit) -> Option<BackgroundModel2D> {
    if let Some(p) = &config.background_model_path {
        // Auto-substitute `lxp1` -> target unit prefix if user passed Lxp10 template path
        let s = p.to_string_lossy();
        let candidate = if s.contains("lxp10") && unit != LaxpcUnit::Lxp10 {
            PathBuf::from(s.replace("lxp10", &format!("{}0", unit.prefix())))
        } else {
            p.clone()
        };
        if candidate.is_file() {
            return BackgroundModel2D::from_file(&candidate).ok();
        }
    }
    None
}

fn accumulate_diag(dst: &mut FrameDiagnostics, src: FrameDiagnostics) {
    dst.total_raw_frames += src.total_raw_frames;
    dst.retained_frames += src.retained_frames;
    dst.rejected_header += src.rejected_header;
    dst.rejected_stray += src.rejected_stray;
    dst.duplicates_replaced += src.duplicates_replaced;
    dst.duplicates_ignored += src.duplicates_ignored;
    dst.out_of_order_count += src.out_of_order_count;
    dst.total_missing_frames += src.total_missing_frames;
    dst.single_frame_gaps += src.single_frame_gaps;
    dst.multi_frame_gaps += src.multi_frame_gaps;
    dst.gaps.extend(src.gaps);
}

/// Build a 1-second cadence `AstroSat.orb` file from the merged MKF (`POSX..VELZ`, matching `laxpcl1.f:orbit`).
pub fn write_orbit_from_mkf(path: &Path, mkf: &MkfTable, meta: &OgipMeta) -> Result<()> {
    let mkf_t = mkf.times();
    if mkf_t.len() < 2 {
        return Err(LaxpcError::Gti("MKF has < 2 rows for orbit generation".into()));
    }
    let px = mkf
        .column_any(&["POSX", "X"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing POSX column".into()))?;
    let py = mkf
        .column_any(&["POSY", "Y"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing POSY column".into()))?;
    let pz = mkf
        .column_any(&["POSZ", "Z"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing POSZ column".into()))?;
    let vx = mkf
        .column_any(&["VELX", "VX"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing VELX column".into()))?;
    let vy = mkf
        .column_any(&["VELY", "VY"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing VELY column".into()))?;
    let vz = mkf
        .column_any(&["VELZ", "VZ"])
        .ok_or_else(|| LaxpcError::Gti("MKF missing VELZ column".into()))?;

    const MAX_ORBIT_SPAN_SEC: f64 = 40_000_000.0; // > 1.25 years at 1 Hz
    let t_first = mkf_t[0];
    let t_last = mkf_t[mkf_t.len() - 1];
    if !t_first.is_finite() || !t_last.is_finite() || t_last < t_first {
        return Err(LaxpcError::Gti(format!(
            "invalid MKF time range for orbit generation: {t_first}..{t_last}"
        )));
    }
    let t_start = t_first.floor() - 3.0;
    let t_stop = t_last.ceil() + 5.0;
    let span = t_stop - t_start;
    if !span.is_finite() || span < 0.0 || span > MAX_ORBIT_SPAN_SEC {
        return Err(LaxpcError::Gti(format!(
            "MKF orbit time span out of bounds ({span} s; max {MAX_ORBIT_SPAN_SEC} s)"
        )));
    }
    let n = (span as usize) + 1;

    let mut times = Vec::with_capacity(n);
    let mut pos_km = Vec::with_capacity(n);
    let mut vel_kms = Vec::with_capacity(n);

    for k in 0..n {
        let t = t_start + (k as f64);
        times.push(t);
        pos_km.push([
            crate::bb::interp_linear(mkf_t, px, t),
            crate::bb::interp_linear(mkf_t, py, t),
            crate::bb::interp_linear(mkf_t, pz, t),
        ]);
        vel_kms.push([
            crate::bb::interp_linear(mkf_t, vx, t),
            crate::bb::interp_linear(mkf_t, vy, t),
            crate::bb::interp_linear(mkf_t, vz, t),
        ]);
    }

    let mut orb_meta = meta.clone();
    orb_meta.tstart = t_start;
    orb_meta.tstop = t_stop;
    orb_meta.detnam = None;
    write_orbit_fits(path, &times, &pos_km, &vel_kms, &orb_meta)
}

/// Run the LAXPC pipeline for a single unit.
pub fn run_unit_pipeline(
    unit: LaxpcUnit,
    orbits: &[OrbitInputFiles],
    config: &PipelineConfig,
) -> Result<UnitPipelineResult> {
    let (res, _ev, mkf, meta) = run_unit_internal(unit, orbits, config)?;
    if config.write_orbit_file {
        let orb_path = config.output_dir.join("AstroSat.orb");
        write_orbit_from_mkf(&orb_path, &mkf, &meta)?;
    }
    Ok(res)
}

/// Run the full multi-unit LAXPC Level-1 -> Level-2 pipeline in parallel using `rayon`.
///
/// Processes all requested units (`LAXPC10`, `LAXPC20`, `LAXPC30`) concurrently, writes
/// Format B per-unit products, `AstroSat.orb`, and (if enabled) Format A combined
/// `level2.event.fits` (`TIME`, `Channel`, `Layer`, `LAXPC_No`).
pub fn run_pipeline(config: &PipelineConfig) -> Result<MultiUnitPipelineResult> {
    if config.units.is_empty() {
        return Err(LaxpcError::Config("no LAXPC units selected".into()));
    }
    fs::create_dir_all(&config.output_dir)?;

    // Resolve orbit file lists for each requested unit
    let mut unit_orbits: Vec<(LaxpcUnit, Vec<OrbitInputFiles>)> = Vec::new();
    for &unit in &config.units {
        let orbits = if let Some(explicit) = config.orbits_by_unit.get(&unit) {
            explicit.clone()
        } else if let Some(root) = &config.input_dir {
            discover_l1_orbits(root, unit)?
        } else {
            return Err(LaxpcError::Config(format!(
                "neither input_dir nor explicit orbits provided for {unit}"
            )));
        };
        unit_orbits.push((unit, orbits));
    }

    // Process all units in parallel with rayon
    let per_unit_res: Result<Vec<InternalUnitRun>> = unit_orbits
        .par_iter()
        .map(|(unit, orbits)| run_unit_internal(*unit, orbits, config))
        .collect();
    let per_unit = per_unit_res?;

    let mut unit_results = Vec::with_capacity(per_unit.len());
    let mut combined_events: Vec<CalibratedEvent> = Vec::new();
    let mut combined_gti_utc = GtiSet::default();
    let mut first_mkf: Option<(MkfTable, OgipMeta)> = None;

    for (res, evs, mkf, meta) in per_unit {
        if first_mkf.is_none() {
            first_mkf = Some((mkf, meta));
        }
        if config.write_format_a_combined {
            combined_gti_utc = if combined_gti_utc.is_empty() {
                res.gti.gti_utc.clone()
            } else {
                combined_gti_utc.union(&res.gti.gti_utc)
            };
            combined_events.extend(evs);
        }
        unit_results.push(res);
    }

    // Write shared AstroSat.orb
    let orbit_fits_path = if config.write_orbit_file {
        if let Some((mkf, meta)) = &first_mkf {
            let p = config.output_dir.join("AstroSat.orb");
            write_orbit_from_mkf(&p, mkf, meta)?;
            Some(p)
        } else {
            None
        }
    } else {
        None
    };

    // Optional Format A combined multi-unit event FITS (`level2.event.fits`)
    let (format_a_event_fits_path, format_a_total_events) = if config.write_format_a_combined {
        // Deterministic multi-core sort with full tie-breakers (`--threads 1` == `--threads N`)
        combined_events.par_sort_by(|a, b| {
            a.time_utc
                .total_cmp(&b.time_utc)
                .then_with(|| a.detector_id.cmp(&b.detector_id))
                .then_with(|| a.anode_id.cmp(&b.anode_id))
                .then_with(|| a.channel.cmp(&b.channel))
                .then_with(|| a.layer.cmp(&b.layer))
        });
        let mut meta_a = first_mkf
            .as_ref()
            .map(|(_, m)| m.clone())
            .unwrap_or_default();
        meta_a.detnam = Some("ALL".to_string());
        if let Some(ts) = combined_gti_utc.tstart() {
            meta_a.tstart = ts;
        }
        if let Some(te) = combined_gti_utc.tstop() {
            meta_a.tstop = te;
        }
        meta_a.exposure = combined_gti_utc.exposure();

        let gti_pairs: Vec<(f64, f64)> = combined_gti_utc
            .intervals()
            .iter()
            .map(|iv| (iv.start, iv.stop))
            .collect();
        let fa_channels = config.spectrum_channels.unwrap_or(1024);
        let path_a = config.output_dir.join("level2.event.fits");
        write_format_a_event_fits_with_channels(
            &path_a,
            &combined_events,
            &gti_pairs,
            fa_channels,
            &meta_a,
        )?;
        (Some(path_a), combined_events.len())
    } else {
        (None, 0)
    };

    Ok(MultiUnitPipelineResult {
        unit_results,
        orbit_fits_path,
        format_a_event_fits_path,
        format_a_total_events,
    })
}
