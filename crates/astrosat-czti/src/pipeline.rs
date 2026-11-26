//! End-to-end AstroSat CZTI Level-1 -> Level-2 pipeline orchestrator (`run_pipeline`)
//! and file-level stage entrypoints for `astrosat-cli`.

use std::fs;
use std::path::{Path, PathBuf};

use astrosat_core::{CztiTct, MkfTable, CZTI_DEFAULT_FEB_LAG};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::bindata::{
    apply_additional_lld, bindata_quad, calculate_renormalized_weights, combine_lightcurves,
    combine_responses, combine_spectra, compute_5x5_avg_open_frac, generate_quad_response,
    radec_to_czti_tx_ty, raytrace_quad_shadow, raytrace_quad_shadow_simple, validate_fov,
    BinDataConfig, BinDataOutputType, CztiLightCurve, CztiResponse, CztiSpectrum,
    DETWISE_TX_SHIFT, DETWISE_TY_SHIFT,
};
use crate::bunchclean::{bunchclean_quad, BunchCleanConfig};
use crate::datasel::datasel_quad;
use crate::dpigen::{generate_dph_dpi, DpiGenConfig};
use crate::error::{CztiError, Result};
use crate::evtclean::{evtclean_quad, EvtCleanConfig};
use crate::evtencor::{evtencor_quad, EvtEnCorConfig};
use crate::fits_io::{
    extract_mkf_rpy_rad, read_aux_tables, read_badpix_fits, read_czti_gtis, read_l1_quad_packets,
    read_l2_event_file, read_mkf_bindata_attitude, write_aspect_file, write_badpix_file,
    write_bunch_file, write_dph_file, write_dpi_file, write_gti_file, write_hdr_file,
    write_l2_event_file, write_lightcurve_file, write_livetime_file, write_response_file,
    write_sky_image_file, write_spectrum_file, CztiCaldb, CztiObsMeta, EventFileWriteOpts,
};
use crate::gaas::compute_quad_aspect;
use crate::gtigen::{
    generate_czti_gtis, intervals_from_mkf_mask, parse_mkf_threshold_text, parse_user_gti_text,
    subtract_intervals, GtiGenConfig,
};
use crate::image::{
    balanced_cross_correlation, balanced_cross_correlation_from_mask, combine_sky_images,
    fft_cross_correlation, ImageAlgorithm, ImageConfig,
};
use crate::noiseclean::{
    eventsep_quad, flickpixclean_quad_with_gti, heavybunchclean_quad_with_exposure,
    noisypixclean_quad, superbunchclean_quad, NoiseCleanConfig,
};
use crate::pha2energy::{pha2energy_quad, QuadTempSeries};
use crate::pixclean::{pixclean_quad, PixCleanConfig};
use crate::science2event::{decode_quad_packets, Science2EventConfig};
use crate::stages::{read_aspect_file, read_bunch_file, read_livetime_file};
use crate::types::{
    detid_pixid_to_detxy, BadpixTable, BunchTable, EventSelectionAlgo, FrameHeaderTable, GtiType,
    QuadAspect, QuadEvents, QuadLivetime, SsmDataTable, TempTable, VetoSpectrumTable, NUM_QUADS,
    PIXELS_PER_QUAD, QUAD_DIM,
};

/// Top-level configuration for running the full AstroSat CZTI Level-1 -> Level-2 pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CztiPipelineConfig {
    /// Path to Level-1 `modeM0` science FITS file (`*cztM0_level1.fits`).
    pub l1_m0_file: PathBuf,
    /// Optional path to Level-1 `modeSS` science FITS file (`*cztSS_level1.fits`) for `TEMP` data.
    pub l1_ss_file: Option<PathBuf>,
    /// Path to Level-1 Time Correlation Table (`*_level1.tct`).
    pub tct_file: PathBuf,
    /// Path to Level-2/Level-1 MKF housekeeping file (`*.mkf`).
    pub mkf_file: PathBuf,
    /// Root path to CZTI CALDB directory (containing `data/as1/czti/caldb.indx` or `caldb.indx`).
    pub caldb_dir: PathBuf,
    /// Output directory where all Level-2 products are written.
    pub output_dir: PathBuf,
    /// Output filename prefix (auto-derived from `l1_m0_file` stem if `None`).
    pub output_prefix: Option<String>,
    /// Optional path to `mkfThresholds.txt` (uses standard CZTI v3.0.1 defaults if `None`).
    pub mkf_threshold_file: Option<PathBuf>,
    /// Optional path to user ASCII GTI file.
    pub user_gti_file: Option<PathBuf>,
    /// Entry stage (`1`, `2`, or `3`; default: `1`).
    pub entry_stage: u8,
    /// Exit stage (`1`, `2`, or `3`; default: `3`).
    pub exit_stage: u8,
    /// GTI mode (`Quad` or `Common`; default: `Quad`).
    pub gti_type: GtiType,
    /// Event selection & cleaning algorithm (`Default` or `Generalized`).
    pub evt_sel_algo: EventSelectionAlgo,
    /// Per-stage tuning configs.
    pub science2event: Science2EventConfig,
    pub bunchclean: BunchCleanConfig,
    pub gtigen: GtiGenConfig,
    pub pixclean: PixCleanConfig,
    pub noiseclean: NoiseCleanConfig,
    pub evtencor: EvtEnCorConfig,
    pub evtclean: EvtCleanConfig,
    pub dpigen: DpiGenConfig,
    pub image: ImageConfig,
    pub bindata: BinDataConfig,
}

impl CztiPipelineConfig {
    /// Convenience constructor with standard CZTI Pipeline v3.0.1 defaults.
    pub fn new(
        l1_m0_file: impl Into<PathBuf>,
        l1_ss_file: Option<PathBuf>,
        tct_file: impl Into<PathBuf>,
        mkf_file: impl Into<PathBuf>,
        caldb_dir: impl Into<PathBuf>,
        output_dir: impl Into<PathBuf>,
    ) -> Self {
        Self {
            l1_m0_file: l1_m0_file.into(),
            l1_ss_file,
            tct_file: tct_file.into(),
            mkf_file: mkf_file.into(),
            caldb_dir: caldb_dir.into(),
            output_dir: output_dir.into(),
            output_prefix: None,
            mkf_threshold_file: None,
            user_gti_file: None,
            entry_stage: 1,
            exit_stage: 3,
            gti_type: GtiType::Quad,
            evt_sel_algo: EventSelectionAlgo::Default,
            science2event: Science2EventConfig::default(),
            bunchclean: BunchCleanConfig::default(),
            gtigen: GtiGenConfig::default(),
            pixclean: PixCleanConfig::default(),
            noiseclean: NoiseCleanConfig::default(),
            evtencor: EvtEnCorConfig::default(),
            evtclean: EvtCleanConfig::default(),
            dpigen: DpiGenConfig::default(),
            image: ImageConfig::default(),
            bindata: BinDataConfig::default(),
        }
    }
}

/// Summary of generated files and per-stage event statistics from `run_pipeline`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CztiPipelineSummary {
    pub obs_id: Option<String>,
    pub stage1_raw_events: [usize; NUM_QUADS],
    pub stage1_bunchclean_events: [usize; NUM_QUADS],
    pub stage2_clean_events: [usize; NUM_QUADS],
    pub quad_exposures_sec: [f64; NUM_QUADS],
    pub quad_aspect_deg: [(f64, f64, f64); NUM_QUADS],
    pub bc_evt_path: Option<PathBuf>,
    pub clean_evt_path: Option<PathBuf>,
    pub gti_path: Option<PathBuf>,
    pub badpix_path: Option<PathBuf>,
    pub livetime_path: Option<PathBuf>,
    pub aspect_paths: Vec<PathBuf>,
    pub dph_path: Option<PathBuf>,
    pub dpi_path: Option<PathBuf>,
    pub image_paths: Vec<PathBuf>,
    pub spectrum_paths: Vec<PathBuf>,
    pub lightcurve_paths: Vec<PathBuf>,
    pub response_paths: Vec<PathBuf>,
    pub weighted_evt_path: Option<PathBuf>,
}

/// Execute the full AstroSat CZTI Level-1 -> Level-2 pipeline in memory with `rayon` parallelism
/// across quadrants, writing all standard CZTI Pipeline v3.0.1 FITS products to `config.output_dir`.
pub fn run_pipeline(config: &CztiPipelineConfig) -> Result<CztiPipelineSummary> {
    if config.entry_stage < 1 || config.exit_stage > 3 || config.entry_stage > config.exit_stage {
        return Err(CztiError::InvalidConfig(format!(
            "Invalid stage range {}..={}",
            config.entry_stage, config.exit_stage
        )));
    }

    fs::create_dir_all(&config.output_dir)?;

    let prefix = config.output_prefix.clone().unwrap_or_else(|| {
        let stem = config
            .l1_m0_file
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "czti".to_string());
        stem.replace("_level1", "_level2")
    });

    let caldb = CztiCaldb::open(&config.caldb_dir)?;
    let gti_word = match config.gti_type {
        GtiType::Quad => "quad",
        GtiType::Common => "common",
    };

    // Determine initial metadata from L1 (entry_stage==1) or existing Stage 1/2 products
    let bc_evt_default = config.output_dir.join(format!("{prefix}_bc.evt"));
    let clean_evt_default = config
        .output_dir
        .join(format!("{prefix}_{gti_word}_clean.evt"));
    let mut meta = if config.entry_stage == 1 {
        CztiObsMeta::from_fits(&config.l1_m0_file)?
    } else if config.entry_stage == 2 && bc_evt_default.exists() {
        CztiObsMeta::from_fits(&bc_evt_default)?
    } else if config.entry_stage == 3 && clean_evt_default.exists() {
        CztiObsMeta::from_fits(&clean_evt_default)?
    } else {
        CztiObsMeta::from_fits(&config.l1_m0_file)?
    };
    meta.caldb_ver = caldb.version().to_string();

    let mut summary = CztiPipelineSummary {
        obs_id: meta.obs_id.clone(),
        ..Default::default()
    };

    let ebounds = caldb.read_ebounds(meta.tstart)?;

    // =========================================================================
    // STAGE 1: science2event (M0 + SS) -> bunchclean -> pha2energy
    // (Skipped when `entry_stage >= 2` by loading existing Stage 1 outputs)
    // =========================================================================
    let (
        mut combined_veto,
        mut combined_ssm,
        mut combined_temp,
        m0_gtis,
        m0_bti,
        bc_events,
        bc_bunches,
        bc_livetimes,
        bc_pix_exp,
        ncompat,
    ) = if config.entry_stage == 1 {
        let tct = CztiTct::from_fits(&config.tct_file)?;
        let (m0_gtis, m0_bti) = read_czti_gtis(&config.l1_m0_file)?;

        // Parallel decode across all 4 quadrants
        let decoded_m0: Vec<_> = (0..NUM_QUADS)
            .into_par_iter()
            .map(|q| -> Result<_> {
                let raw = read_l1_quad_packets(&config.l1_m0_file, q)?;
                let stream = decode_quad_packets(q as u8, &raw, &config.science2event, |czt_sec| {
                    tct.interpolate(czt_sec, CZTI_DEFAULT_FEB_LAG, None)
                        .unwrap_or(czt_sec)
                });
                Ok(stream)
            })
            .collect::<Result<Vec<_>>>()?;

        // Optionally decode modeSS for TEMP & SSM Data
        let mut combined_temp = TempTable::default();
        let mut combined_ssm = SsmDataTable::default();
        if let Some(ref ss_path) = config.l1_ss_file {
            if ss_path.exists() {
                let decoded_ss: Vec<_> = (0..NUM_QUADS)
                    .into_par_iter()
                    .map(|q| -> Result<_> {
                        let raw = read_l1_quad_packets(ss_path, q)?;
                        Ok(decode_quad_packets(
                            q as u8,
                            &raw,
                            &config.science2event,
                            |czt_sec| {
                                tct.interpolate(czt_sec, CZTI_DEFAULT_FEB_LAG, None)
                                    .unwrap_or(czt_sec)
                            },
                        ))
                    })
                    .collect::<Result<Vec<_>>>()?;
                for s in decoded_ss {
                    combined_temp.time.extend(s.temp.time);
                    combined_temp.cztseccnt.extend(s.temp.cztseccnt);
                    combined_temp.temperature.extend(s.temp.temperature);
                    combined_temp.quad_flag.extend(s.temp.quad_flag);

                    combined_ssm.time.extend(s.ssm_data.time);
                    combined_ssm.cztseccnt.extend(s.ssm_data.cztseccnt);
                    combined_ssm.veto_spec.extend(s.ssm_data.veto_spec);
                    combined_ssm.czt_spec.extend(s.ssm_data.czt_spec);
                    combined_ssm.czt_spec_veto.extend(s.ssm_data.czt_spec_veto);
                    combined_ssm.czt_spec_alpha.extend(s.ssm_data.czt_spec_alpha);
                    combined_ssm.quad_id.extend(s.ssm_data.quad_id);
                }
            }
        }

        let mut combined_veto = VetoSpectrumTable::default();
        let mut combined_hdr = FrameHeaderTable::default();
        for (q, d) in decoded_m0.iter().enumerate() {
            summary.stage1_raw_events[q] = d.events.len();
            combined_veto.time.extend_from_slice(&d.veto_spec.time);
            combined_veto.cztseccnt.extend_from_slice(&d.veto_spec.cztseccnt);
            combined_veto.veto_spec.extend_from_slice(&d.veto_spec.veto_spec);
            combined_veto.quad_id.extend_from_slice(&d.veto_spec.quad_id);

            combined_hdr.time.extend_from_slice(&d.headers.time);
            combined_hdr.cztseccnt.extend_from_slice(&d.headers.cztseccnt);
            combined_hdr.data_id.extend_from_slice(&d.headers.data_id);
            combined_hdr.mode_id.extend_from_slice(&d.headers.mode_id);
            combined_hdr.word_count.extend_from_slice(&d.headers.word_count);
            combined_hdr.frame_no.extend_from_slice(&d.headers.frame_no);
            combined_hdr.status_flag.extend_from_slice(&d.headers.status_flag);
            combined_hdr.write_pkt_no.extend_from_slice(&d.headers.write_pkt_no);
            combined_hdr.read_pkt_no.extend_from_slice(&d.headers.read_pkt_no);
            combined_hdr.command.extend_from_slice(&d.headers.command);
            combined_hdr.laxpc_time.extend_from_slice(&d.headers.laxpc_time);
            combined_hdr.dcnt.extend_from_slice(&d.headers.dcnt);
            combined_hdr.error_count.extend_from_slice(&d.headers.error_count);
            combined_hdr.error_flag.extend_from_slice(&d.headers.error_flag);
            combined_hdr.boot_page_no.extend_from_slice(&d.headers.boot_page_no);
            combined_hdr.czt_no.extend_from_slice(&d.headers.czt_no);
            combined_hdr.czt_status.extend_from_slice(&d.headers.czt_status);
            combined_hdr.evt_read_mode.extend_from_slice(&d.headers.evt_read_mode);
            combined_hdr.cmd_status.extend_from_slice(&d.headers.cmd_status);
            combined_hdr.buffer_no.extend_from_slice(&d.headers.buffer_no);
            combined_hdr.base_add.extend_from_slice(&d.headers.base_add);
            combined_hdr.veto_spec_range.extend_from_slice(&d.headers.veto_spec_range);
            combined_hdr.channel_no.extend_from_slice(&d.headers.channel_no);
            combined_hdr.adc_output.extend_from_slice(&d.headers.adc_output);
            combined_hdr.cmd_1sec.extend_from_slice(&d.headers.cmd_1sec);
            combined_hdr.alpha_count.extend_from_slice(&d.headers.alpha_count);
            combined_hdr.veto_count.extend_from_slice(&d.headers.veto_count);
            combined_hdr.czt_count_lt_uld.extend_from_slice(&d.headers.czt_count_lt_uld);
            combined_hdr.czt_count_gt_uld.extend_from_slice(&d.headers.czt_count_gt_uld);
            combined_hdr.czt_data_read.extend_from_slice(&d.headers.czt_data_read);
        }

        // Populate TSTART/TSTOP if absent in primary header
        if meta.tstart == 0.0 {
            meta.tstart = decoded_m0
                .iter()
                .filter_map(|d| d.events.time.first().copied())
                .fold(f64::INFINITY, f64::min);
            if !meta.tstart.is_finite() {
                meta.tstart = 0.0;
            }
        }
        if meta.tstop == 0.0 {
            meta.tstop = decoded_m0
                .iter()
                .filter_map(|d| d.events.time.last().copied())
                .fold(0.0f64, f64::max);
        }

        let hdr_path = config.output_dir.join(format!("{prefix}.hdr"));
        write_hdr_file(&hdr_path, &meta, &combined_hdr)?;

        // Bunch cleaning in parallel across quadrants
        let bc_cfg = match config.evt_sel_algo {
            EventSelectionAlgo::Generalized => BunchCleanConfig::generalized_compat(),
            EventSelectionAlgo::Default => config.bunchclean.clone(),
        };
        let bc_results: Vec<_> = (0..NUM_QUADS)
            .into_par_iter()
            .map(|q| {
                bunchclean_quad(
                    &decoded_m0[q].events,
                    &decoded_m0[q].bunches,
                    &m0_gtis.quad[q],
                    &bc_cfg,
                )
            })
            .collect();

        let ncompat = bc_results[0].ncompat;
        let mut bc_events: [QuadEvents; NUM_QUADS] = [
            bc_results[0].events.clone(),
            bc_results[1].events.clone(),
            bc_results[2].events.clone(),
            bc_results[3].events.clone(),
        ];
        let bc_bunches: [BunchTable; NUM_QUADS] = [
            bc_results[0].out_bunches.clone(),
            bc_results[1].out_bunches.clone(),
            bc_results[2].out_bunches.clone(),
            bc_results[3].out_bunches.clone(),
        ];
        let bc_livetimes: [QuadLivetime; NUM_QUADS] = [
            bc_results[0].livetime.clone(),
            bc_results[1].livetime.clone(),
            bc_results[2].livetime.clone(),
            bc_results[3].livetime.clone(),
        ];
        let bc_pix_exp: [Vec<f64>; NUM_QUADS] = [
            bc_results[0].pixel_exposure.clone(),
            bc_results[1].pixel_exposure.clone(),
            bc_results[2].pixel_exposure.clone(),
            bc_results[3].pixel_exposure.clone(),
        ];
        for q in 0..NUM_QUADS {
            summary.stage1_bunchclean_events[q] = bc_events[q].len();
        }

        let bunch_path = config.output_dir.join(format!("{prefix}_bc_bunch.fits"));
        write_bunch_file(&bunch_path, &meta, &bc_bunches)?;
        let bc_lv_path = config.output_dir.join(format!("{prefix}_bc_livetime.fits"));
        write_livetime_file(&bc_lv_path, &meta, &bc_livetimes)?;

        // PHA -> Energy & PI via CALDB GAIN + EBOUNDS (parallelized across quadrants)
        let gain_offsets: Vec<_> = (0..NUM_QUADS)
            .map(|q| caldb.read_gain_offset(q, meta.tstart))
            .collect::<Result<Vec<_>>>()?;
        let temp_series: Vec<QuadTempSeries> = (0..NUM_QUADS)
            .map(|q| QuadTempSeries::from_temp_table(&combined_temp, q as u8))
            .collect();

        bc_events
            .par_iter_mut()
            .enumerate()
            .for_each(|(q, ev)| pha2energy_quad(ev, &temp_series[q], &gain_offsets[q], &ebounds));

        write_l2_event_file(
            &bc_evt_default,
            &meta,
            &bc_events,
            &combined_veto,
            &combined_ssm,
            &combined_temp,
            &m0_gtis,
            &EventFileWriteOpts {
                gti_type: None,
                ncompat: Some(ncompat),
                gaincor_applied: false,
                lld_cut: false,
                uld_cut: false,
                pixel_exposure: Some(bc_pix_exp.clone()),
            },
        )?;
        summary.bc_evt_path = Some(bc_evt_default.clone());
        (
            combined_veto,
            combined_ssm,
            combined_temp,
            m0_gtis,
            m0_bti,
            bc_events,
            bc_bunches,
            bc_livetimes,
            bc_pix_exp,
            ncompat,
        )
    } else if config.entry_stage == 2 {
        if !bc_evt_default.exists() {
            return Err(CztiError::MissingFile(bc_evt_default));
        }
        let (_, bc_events, m0_gtis, pexp_opt, _) = read_l2_event_file(&bc_evt_default)?;
        let (_, m0_bti) = read_czti_gtis(&bc_evt_default)?;
        let (v, s, t) = read_aux_tables(&bc_evt_default);
        let bunch_path = config.output_dir.join(format!("{prefix}_bc_bunch.fits"));
        let bc_bunches = if bunch_path.exists() {
            read_bunch_file(&bunch_path)?
        } else {
            [
                BunchTable::new(0),
                BunchTable::new(1),
                BunchTable::new(2),
                BunchTable::new(3),
            ]
        };
        let bc_lv_path = config.output_dir.join(format!("{prefix}_bc_livetime.fits"));
        let bc_livetimes = if bc_lv_path.exists() {
            read_livetime_file(&bc_lv_path)?
        } else {
            [
                QuadLivetime::default(),
                QuadLivetime::default(),
                QuadLivetime::default(),
                QuadLivetime::default(),
            ]
        };
        let bc_pix_exp = pexp_opt.unwrap_or_else(|| {
            [
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
            ]
        });
        for q in 0..NUM_QUADS {
            summary.stage1_bunchclean_events[q] = bc_events[q].len();
        }
        summary.bc_evt_path = Some(bc_evt_default.clone());
        let ncompat = if config.evt_sel_algo == EventSelectionAlgo::Generalized {
            1
        } else {
            0
        };
        (
            v,
            s,
            t,
            m0_gtis,
            m0_bti,
            bc_events,
            bc_bunches,
            bc_livetimes,
            bc_pix_exp,
            ncompat,
        )
    } else {
        // entry_stage == 3: dummy placeholders; Stage 2 block loads existing clean products
        (
            VetoSpectrumTable::default(),
            SsmDataTable::default(),
            TempTable::default(),
            crate::gtigen::CztiGtiOutput::default(),
            Vec::new(),
            [
                QuadEvents::new(0),
                QuadEvents::new(1),
                QuadEvents::new(2),
                QuadEvents::new(3),
            ],
            [
                BunchTable::new(0),
                BunchTable::new(1),
                BunchTable::new(2),
                BunchTable::new(3),
            ],
            [
                QuadLivetime::default(),
                QuadLivetime::default(),
                QuadLivetime::default(),
                QuadLivetime::default(),
            ],
            [
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
            ],
            0,
        )
    };

    if config.exit_stage < 2 {
        return Ok(summary);
    }

    // =========================================================================
    // STAGE 2: gtigen -> gaas -> datasel -> pixclean/noiseclean -> evtencor -> evtclean
    // (Skipped when `entry_stage == 3` by loading existing Stage 2 outputs)
    // =========================================================================
    let gaincor_caldb: Vec<_> = (0..NUM_QUADS)
        .map(|q| caldb.read_gaincor_lld_uld(q, meta.tstart))
        .collect::<Result<Vec<_>>>()?;

    let (
        final_gtis,
        aspects,
        final_badpix,
        final_pix_exp,
        clean_quads_arr,
        final_lv_arr,
    ) = if config.entry_stage <= 2 {
        let mkf = MkfTable::from_fits(&config.mkf_file)?;
        let mkf_rules = if let Some(ref tf) = config.mkf_threshold_file {
            let txt = fs::read_to_string(tf)?;
            parse_mkf_threshold_text(&txt)
        } else {
            config.gtigen.thresholds.clone()
        };

        let n_mkf = mkf.len();
        let mut mkf_good = vec![true; n_mkf];
        for rule in &mkf_rules {
            if !rule.enabled {
                continue;
            }
            if let Some(col) = mkf.column(&rule.column) {
                for (i, &val) in col.iter().enumerate() {
                    if !val.is_finite() {
                        mkf_good[i] = false;
                        continue;
                    }
                    if let Some(lo) = rule.min {
                        if val < lo {
                            mkf_good[i] = false;
                        }
                    }
                    if let Some(hi) = rule.max {
                        if val > hi {
                            mkf_good[i] = false;
                        }
                    }
                }
            }
        }
        let mkf_gti = intervals_from_mkf_mask(mkf.times(), &mkf_good);

        let user_gti = if let Some(ref ug_path) = config.user_gti_file {
            let txt = fs::read_to_string(ug_path)?;
            Some(parse_user_gti_text(&txt)?)
        } else {
            None
        };

        let mut final_gtis = generate_czti_gtis(
            &mkf_gti,
            &m0_gtis.common,
            &m0_gtis.quad,
            if m0_bti.is_empty() { None } else { Some(&m0_bti) },
            user_gti.as_deref(),
            config.gtigen.min_exposure,
        )?;

        // GAAS: compute per-quadrant average aspect (RA, DEC, TWIST)
        let mkf_rpy = extract_mkf_rpy_rad(&mkf)?;
        let aspects: Vec<QuadAspect> = (0..NUM_QUADS)
            .map(|q| -> Result<QuadAspect> {
                let (tname, tmat) = caldb.read_teldef(q, meta.tstart)?;
                let qgti = match config.gti_type {
                    GtiType::Common => &final_gtis.common,
                    GtiType::Quad => &final_gtis.quad[q],
                };
                let t0 = if bc_events[q].tstart > 0.0 {
                    bc_events[q].tstart
                } else {
                    meta.tstart
                };
                let t1 = if bc_events[q].tstop > 0.0 {
                    bc_events[q].tstop
                } else {
                    meta.tstop
                };
                Ok(compute_quad_aspect(
                    q as u8,
                    tname,
                    &tmat,
                    &mkf_rpy,
                    t0,
                    t1,
                    if qgti.is_empty() {
                        None
                    } else {
                        Some(qgti.as_slice())
                    },
                ))
            })
            .collect::<Result<Vec<_>>>()?;

        for (q, asp) in aspects.iter().enumerate() {
            summary.quad_aspect_deg[q] = (asp.ra_deg, asp.dec_deg, asp.twist_deg);
            let apath = config.output_dir.join(format!("{prefix}.aspect_Q{q}"));
            write_aspect_file(&apath, &meta, asp)?;
            summary.aspect_paths.push(apath);
        }

        // Load CALDB badpix
        let caldb_badpix = caldb.read_badpix(meta.tstart)?;

        let mut final_badpix = BadpixTable::default();
        let mut final_livetimes: Vec<QuadLivetime> = Vec::with_capacity(NUM_QUADS);
        let mut final_pix_exp: [Vec<f64>; NUM_QUADS] = [
            vec![0.0f64; PIXELS_PER_QUAD],
            vec![0.0f64; PIXELS_PER_QUAD],
            vec![0.0f64; PIXELS_PER_QUAD],
            vec![0.0f64; PIXELS_PER_QUAD],
        ];
        let mut clean_quads: Vec<QuadEvents> = Vec::with_capacity(NUM_QUADS);
        let mut double_quads: Vec<QuadEvents> = Vec::with_capacity(NUM_QUADS);

        let write_dbl = match config.evt_sel_algo {
            EventSelectionAlgo::Default => config.pixclean.write_double_evt,
            EventSelectionAlgo::Generalized => config.noiseclean.write_double_evt,
        };
        let dbl_clean_cfg = EvtCleanConfig {
            is_double_event: true,
            ..config.evtclean.clone()
        };

        let stage2_per_quad: Vec<_> = (0..NUM_QUADS)
            .into_par_iter()
            .map(|q| {
                let qgti = match config.gti_type {
                    GtiType::Common => &final_gtis.common,
                    GtiType::Quad => &final_gtis.quad[q],
                };
                let (mut single_ev, mut dbl_ev, lv, bp_flags, pexp, sb_btis) =
                    match config.evt_sel_algo {
                        EventSelectionAlgo::Default => {
                            let ds = datasel_quad(&bc_events[q], qgti);
                            let pc = pixclean_quad(
                                &ds,
                                &bc_livetimes[q],
                                qgti,
                                &caldb_badpix.flags[q],
                                &bc_pix_exp[q],
                                &config.pixclean,
                            );
                            (
                                pc.single_events,
                                pc.double_events,
                                pc.livetime,
                                pc.badpix_flags,
                                pc.pixel_exposure,
                                Vec::new(),
                            )
                        }
                        EventSelectionAlgo::Generalized => {
                            let np = noisypixclean_quad(
                                &bc_events[q],
                                &caldb_badpix.flags[q],
                                &bc_pix_exp[q],
                                &config.noiseclean,
                            );
                            let ds = datasel_quad(&np.events, qgti);
                            let sb = superbunchclean_quad(
                                &ds,
                                &bc_bunches[q],
                                &bc_livetimes[q],
                                &config.noiseclean,
                            );
                            let (hb_ev, hb_pexp) = heavybunchclean_quad_with_exposure(
                                &sb.events,
                                &bc_bunches[q],
                                &gaincor_caldb[q].lld,
                                &np.pixel_exposure,
                                &config.noiseclean,
                            );
                            let fp = flickpixclean_quad_with_gti(
                                &hb_ev,
                                &np.badpix_flags,
                                &hb_pexp,
                                Some(qgti),
                                &config.noiseclean,
                            );
                            let (mut singles, mut doubles) =
                                eventsep_quad(&fp.events, config.noiseclean.double_time_sec);
                            if !sb.bti_intervals.is_empty() {
                                let updated_gti = subtract_intervals(qgti, &sb.bti_intervals);
                                let new_exp: f64 =
                                    updated_gti.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
                                singles.exposure = new_exp;
                                doubles.exposure = new_exp;
                            }
                            (
                                singles,
                                doubles,
                                sb.livetime,
                                fp.badpix_flags,
                                fp.pixel_exposure,
                                sb.bti_intervals,
                            )
                        }
                    };

                // evtencor (gain correction + secular drift + LLD/ULD cut)
                let ec = evtencor_quad(
                    &single_ev,
                    &gaincor_caldb[q],
                    meta.tstart,
                    false,
                    &config.evtencor,
                );
                single_ev = ec.events;

                // evtclean (alpha & veto selection)
                let cleaned = evtclean_quad(&single_ev, &config.evtclean);

                // Match legacy `cztpipeline.py:613-637`: also calibrate & clean double events
                if write_dbl && !dbl_ev.is_empty() {
                    let ec_dbl = evtencor_quad(
                        &dbl_ev,
                        &gaincor_caldb[q],
                        meta.tstart,
                        false,
                        &config.evtencor,
                    );
                    dbl_ev = evtclean_quad(&ec_dbl.events, &dbl_clean_cfg);
                }
                (cleaned, dbl_ev, lv, bp_flags, pexp, sb_btis)
            })
            .collect();

        for (q, (cleaned, dbl_ev, lv, bp_flags, pexp, sb_btis)) in
            stage2_per_quad.into_iter().enumerate()
        {
            if !sb_btis.is_empty() {
                final_gtis.quad[q] = subtract_intervals(&final_gtis.quad[q], &sb_btis);
            }
            summary.stage2_clean_events[q] = cleaned.len();
            summary.quad_exposures_sec[q] = cleaned.exposure;
            final_badpix.flags[q] = bp_flags;
            final_pix_exp[q] = pexp;
            final_livetimes.push(lv);
            clean_quads.push(cleaned);
            double_quads.push(dbl_ev);
        }

        let gti_path = config.output_dir.join(format!("{prefix}.gti"));
        write_gti_file(&gti_path, &meta, &final_gtis)?;
        summary.gti_path = Some(gti_path);

        let mut cq_it = clean_quads.into_iter();
        let clean_quads_arr: [QuadEvents; NUM_QUADS] = [
            cq_it.next().unwrap(),
            cq_it.next().unwrap(),
            cq_it.next().unwrap(),
            cq_it.next().unwrap(),
        ];
        let mut lv_it = final_livetimes.into_iter();
        let final_lv_arr: [QuadLivetime; NUM_QUADS] = [
            lv_it.next().unwrap(),
            lv_it.next().unwrap(),
            lv_it.next().unwrap(),
            lv_it.next().unwrap(),
        ];

        let badpix_path = config
            .output_dir
            .join(format!("{prefix}_{gti_word}_badpix.fits"));
        write_badpix_file(&badpix_path, &meta, &final_badpix)?;
        summary.badpix_path = Some(badpix_path);

        let lv_path = config
            .output_dir
            .join(format!("{prefix}_{gti_word}_livetime.fits"));
        write_livetime_file(&lv_path, &meta, &final_lv_arr)?;
        summary.livetime_path = Some(lv_path);

        write_l2_event_file(
            &clean_evt_default,
            &meta,
            &clean_quads_arr,
            &combined_veto,
            &combined_ssm,
            &combined_temp,
            &final_gtis,
            &EventFileWriteOpts {
                gti_type: Some(config.gti_type),
                ncompat: Some(ncompat),
                gaincor_applied: true,
                lld_cut: config.evtencor.apply_lld,
                uld_cut: config.evtencor.apply_uld,
                pixel_exposure: Some(final_pix_exp.clone()),
            },
        )?;
        summary.clean_evt_path = Some(clean_evt_default.clone());

        if write_dbl {
            let mut dq_it = double_quads.into_iter();
            let dbl_arr: [QuadEvents; NUM_QUADS] = [
                dq_it.next().unwrap(),
                dq_it.next().unwrap(),
                dq_it.next().unwrap(),
                dq_it.next().unwrap(),
            ];
            let dbl_path = config
                .output_dir
                .join(format!("{prefix}_{gti_word}_double.evt"));
            write_l2_event_file(
                &dbl_path,
                &meta,
                &dbl_arr,
                &combined_veto,
                &combined_ssm,
                &combined_temp,
                &final_gtis,
                &EventFileWriteOpts {
                    gti_type: Some(config.gti_type),
                    ncompat: Some(ncompat),
                    gaincor_applied: true,
                    lld_cut: config.evtencor.apply_lld,
                    uld_cut: config.evtencor.apply_uld,
                    pixel_exposure: Some(final_pix_exp.clone()),
                },
            )?;
        }

        (
            final_gtis,
            aspects,
            final_badpix,
            final_pix_exp,
            clean_quads_arr,
            final_lv_arr,
        )
    } else {
        // entry_stage == 3: load existing Stage 2 outputs from `output_dir`
        if !clean_evt_default.exists() {
            return Err(CztiError::MissingFile(clean_evt_default));
        }
        let (_, clean_quads_arr, final_gtis, pexp_opt, _) =
            read_l2_event_file(&clean_evt_default)?;
        let (v, s, t) = read_aux_tables(&clean_evt_default);
        combined_veto = v;
        combined_ssm = s;
        combined_temp = t;
        let badpix_path = config
            .output_dir
            .join(format!("{prefix}_{gti_word}_badpix.fits"));
        let final_badpix = if badpix_path.exists() {
            read_badpix_fits(&badpix_path)?
        } else {
            caldb.read_badpix(meta.tstart)?
        };
        let lv_path = config
            .output_dir
            .join(format!("{prefix}_{gti_word}_livetime.fits"));
        let final_lv_arr = if lv_path.exists() {
            read_livetime_file(&lv_path)?
        } else {
            [
                QuadLivetime::default(),
                QuadLivetime::default(),
                QuadLivetime::default(),
                QuadLivetime::default(),
            ]
        };
        let mut aspects = Vec::with_capacity(NUM_QUADS);
        for q in 0..NUM_QUADS {
            let apath = config.output_dir.join(format!("{prefix}.aspect_Q{q}"));
            let asp = if apath.exists() {
                read_aspect_file(&apath, q as u8)?
            } else {
                QuadAspect {
                    quad_id: q as u8,
                    ..Default::default()
                }
            };
            summary.quad_aspect_deg[q] = (asp.ra_deg, asp.dec_deg, asp.twist_deg);
            summary.stage2_clean_events[q] = clean_quads_arr[q].len();
            summary.quad_exposures_sec[q] = clean_quads_arr[q].exposure;
            aspects.push(asp);
        }
        let final_pix_exp = pexp_opt.unwrap_or_else(|| {
            [
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
                vec![1.0; PIXELS_PER_QUAD],
            ]
        });
        summary.clean_evt_path = Some(clean_evt_default);
        (
            final_gtis,
            aspects,
            final_badpix,
            final_pix_exp,
            clean_quads_arr,
            final_lv_arr,
        )
    };

    if config.exit_stage < 3 {
        return Ok(summary);
    }

    // =========================================================================
    // STAGE 3: dpigen -> image -> bindata (.pha, .lc, .rsp, optional weighted .evt)
    // =========================================================================
    let eff_areas: [_; NUM_QUADS] = [
        caldb.read_eff_area(0, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(1, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(2, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(3, meta.tstart).unwrap_or_default(),
    ];
    let dph_dpi = generate_dph_dpi(&clean_quads_arr, &final_badpix, &eff_areas, &config.dpigen);

    let dph_path = config
        .output_dir
        .join(format!("{prefix}_{gti_word}_clean.dph"));
    let dpi_path = config
        .output_dir
        .join(format!("{prefix}_{gti_word}_clean.dpi"));
    write_dph_file(
        &dph_path,
        &meta,
        config.gti_type,
        &dph_dpi.quad_dph,
        &dph_dpi.full_dph,
    )?;
    write_dpi_file(
        &dpi_path,
        &meta,
        config.gti_type,
        &dph_dpi.quad_dpi,
        &dph_dpi.full_dpi,
    )?;
    summary.dph_path = Some(dph_path);
    summary.dpi_path = Some(dpi_path);

    // Sky image reconstruction (`balcc` or `fft`) — load gzipped SHADOWLIB once across quadrants!
    let shadowlib_all = if config.image.algorithm == ImageAlgorithm::BalCC {
        caldb.read_shadowlib_all(meta.tstart).ok()
    } else {
        None
    };
    let mut quad_sky_imgs = Vec::new();
    for q in 0..NUM_QUADS {
        if !config.image.quads.contains(&q) {
            continue;
        }
        let sky_img = match config.image.algorithm {
            ImageAlgorithm::BalCC => {
                if let Some((grid_dim, ref all_shadows)) = shadowlib_all {
                    let refs: Vec<&[f32]> = all_shadows[q].iter().map(|v| v.as_slice()).collect();
                    balanced_cross_correlation(&dph_dpi.quad_dpi[q], &refs, grid_dim, &aspects[q])
                } else {
                    let mask = caldb
                        .read_mask_pattern(q, meta.tstart)
                        .unwrap_or_else(|_| vec![vec![0.5; QUAD_DIM]; QUAD_DIM]);
                    balanced_cross_correlation_from_mask(&dph_dpi.quad_dpi[q], &mask, 15, &aspects[q])
                }
            }
            ImageAlgorithm::Fft => {
                let mask = caldb
                    .read_mask_pattern(q, meta.tstart)
                    .unwrap_or_else(|_| vec![vec![0.5; QUAD_DIM]; QUAD_DIM]);
                fft_cross_correlation(
                    &dph_dpi.quad_dpi[q],
                    &mask,
                    config.image.oversampling_factor,
                    &aspects[q],
                )
            }
        };
        if config.gti_type == GtiType::Quad {
            let img_path = config
                .output_dir
                .join(format!("{prefix}_{gti_word}_clean.img_Q{q}"));
            write_sky_image_file(&img_path, &format!("Q{q}"), &meta, &sky_img)?;
            summary.image_paths.push(img_path);
        }
        quad_sky_imgs.push(sky_img);
    }
    if config.gti_type == GtiType::Common && !quad_sky_imgs.is_empty() {
        let combined = combine_sky_images(&quad_sky_imgs, &aspects[0]);
        let img_path = config
            .output_dir
            .join(format!("{prefix}_{gti_word}_clean.img"));
        write_sky_image_file(&img_path, "IMAGE", &meta, &combined)?;
        summary.image_paths.push(img_path);
    }
    drop(shadowlib_all);

    // Load shared PIXRESP once if response generation is enabled
    let shared_pixresp = if config.bindata.gen_resp
        && matches!(
            config.bindata.output_type,
            BinDataOutputType::Spec | BinDataOutputType::Both
        ) {
        caldb.read_pixresp(meta.tstart).ok()
    } else {
        None
    };

    // Bindata: camera (tx, ty) -> coded-mask ray-tracing -> BKGSCAL weights -> .pha / .lc / .rsp
    let (ra_pnt, dec_pnt, att_t, roll_ra, roll_dec, roll_rot) =
        read_mkf_bindata_attitude(&config.mkf_file)?;
    let (ra_use, dec_use) = config.bindata.resolve_source_radec(ra_pnt, dec_pnt);

    let att_tx_ty: Vec<(f64, f64)> = (0..att_t.len())
        .into_par_iter()
        .map(|i| radec_to_czti_tx_ty(roll_ra[i], roll_dec[i], roll_rot[i], ra_use, dec_use))
        .collect();

    let mut weighted_quads: [QuadEvents; NUM_QUADS] = clean_quads_arr.clone();
    let mut wrote_weighted = false;
    let mut common_specs: Vec<CztiSpectrum> = Vec::new();
    let mut common_lcs: Vec<CztiLightCurve> = Vec::new();
    let mut common_rsps: Vec<CztiResponse> = Vec::new();

    for &q in &config.bindata.quads {
        if q >= NUM_QUADS {
            continue;
        }
        let qgti = match config.gti_type {
            GtiType::Common => &final_gtis.common,
            GtiType::Quad => &final_gtis.quad[q],
        };
        let mut sum_tx = 0.0f64;
        let mut sum_ty = 0.0f64;
        let mut gtisel_tx_ty: Vec<(f64, f64)> = Vec::new();
        for (i, &t) in att_t.iter().enumerate() {
            let in_g = qgti.is_empty() || qgti.iter().any(|&(a, b)| t >= a && t <= b);
            if in_g {
                let (tx, ty) = att_tx_ty[i];
                if config.bindata.mask_weight || config.bindata.gen_resp {
                    validate_fov(tx, ty)?;
                }
                sum_tx += tx;
                sum_ty += ty;
                gtisel_tx_ty.push((tx, ty));
            }
        }
        let (avg_tx, avg_ty) = if !gtisel_tx_ty.is_empty() {
            (
                sum_tx / gtisel_tx_ty.len() as f64,
                sum_ty / gtisel_tx_ty.len() as f64,
            )
        } else if let Some(&(tx, ty)) = att_tx_ty.first() {
            (tx, ty)
        } else {
            (0.0, 0.0)
        };

        // Apply `additional_LLD` per-module floor (`cztbindata.cpp:1930-1931`) for bindata & rspgen
        let mut gc_bind = gaincor_caldb[q].clone();
        apply_additional_lld(q, &mut gc_bind.lld);

        // Build [detx][dety] LLD/ULD maps and pix_exp[detx*64 + dety]
        let mut lld_detmap = [[0i16; QUAD_DIM]; QUAD_DIM];
        let mut uld_detmap = [[511i16; QUAD_DIM]; QUAD_DIM];
        let mut pix_exp_detxy = vec![0.0f64; PIXELS_PER_QUAD];
        for det in 0..16u8 {
            for pix in 0..=255u8 {
                let row = (det as usize) * 256 + (pix as usize);
                let (dx, dy) = detid_pixid_to_detxy(q as u8, det, pix);
                lld_detmap[dx as usize][dy as usize] = gc_bind.lld[row];
                uld_detmap[dx as usize][dy as usize] = gc_bind.uld[row];
                pix_exp_detxy[(dx as usize) * QUAD_DIM + (dy as usize)] = final_pix_exp[q][row];
            }
        }

        // Ray-trace coded-mask shadow & compute BKGSCAL-renormalized mask weights
        let mut tx_sh = DETWISE_TX_SHIFT[q];
        let mut ty_sh = DETWISE_TY_SHIFT[q];
        for m in 0..16 {
            tx_sh[m] += config.bindata.delta_tx;
            ty_sh[m] += config.bindata.delta_ty;
        }
        let mask_ov_opt = caldb.read_mask_oversampled(q, meta.tstart).ok();
        let mask_64_opt = if mask_ov_opt.is_none() {
            Some(
                caldb
                    .read_mask_pattern(q, meta.tstart)
                    .unwrap_or_else(|_| vec![vec![0.5; QUAD_DIM]; QUAD_DIM]),
            )
        } else {
            None
        };
        let (open_frac, raw_w) = if let Some(ref mask_ov) = mask_ov_opt {
            raytrace_quad_shadow(avg_tx, avg_ty, mask_ov, &tx_sh, &ty_sh)
        } else {
            raytrace_quad_shadow_simple(avg_tx, avg_ty, mask_64_opt.as_ref().unwrap())
        };

        let renorm_weights = if config.bindata.mask_weight {
            let bkgscal = caldb
                .read_bkgscal(q, meta.tstart)
                .unwrap_or_else(|_| vec![1.0; PIXELS_PER_QUAD * 512]);
            Some(calculate_renormalized_weights(
                q,
                &open_frac,
                &raw_w,
                &final_badpix,
                config.bindata.badpix_threshold,
                &bkgscal,
                &lld_detmap,
                &uld_detmap,
                &pix_exp_detxy,
                config.bindata.min_active_pix_per_module,
            ))
        } else {
            None
        };

        let syserr = caldb
            .read_syserr(q, meta.tstart)
            .unwrap_or([0.0f32; 512]);

        let rsp_filename = if config.bindata.gen_resp {
            if config.gti_type == GtiType::Common {
                format!("{prefix}_{gti_word}_clean.rsp")
            } else {
                format!("{prefix}_{gti_word}_clean_Q{q}.rsp")
            }
        } else {
            "NONE".to_string()
        };

        let (spec, lc, wevt) = bindata_quad(
            &clean_quads_arr[q],
            &final_badpix,
            &final_lv_arr[q],
            renorm_weights.as_ref(),
            &gc_bind,
            true,
            &syserr,
            meta.tstart,
            meta.tstop,
            avg_tx,
            avg_ty,
            rsp_filename.clone(),
            &config.bindata,
        );

        let quad_rsp = if config.bindata.gen_resp
            && matches!(
                config.bindata.output_type,
                BinDataOutputType::Spec | BinDataOutputType::Both
            ) {
            if let Some(ref pixresp) = shared_pixresp {
                if let Ok(mut rsp_caldb) =
                    caldb.read_rsp_caldb_with_pixresp(q, meta.tstart, pixresp)
                {
                    apply_additional_lld(q, &mut rsp_caldb.lld);
                    // 5x5 `(tx, ty)` attitude-histogram average open fraction (`cztrspgen.cpp:441-505`)
                    let avg_of_5x5 = compute_5x5_avg_open_frac(
                        &gtisel_tx_ty,
                        mask_ov_opt.as_deref(),
                        mask_64_opt.as_deref(),
                        &tx_sh,
                        &ty_sh,
                        &open_frac,
                    );
                    Some(generate_quad_response(
                        q as u8,
                        &rsp_caldb,
                        &final_badpix,
                        config.bindata.badpix_threshold,
                        &avg_of_5x5,
                        renorm_weights.as_ref(),
                        &pix_exp_detxy,
                        ebounds.e_min.clone(),
                        ebounds.e_max.clone(),
                    ))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };

        if config.gti_type == GtiType::Common {
            common_specs.push(spec);
            common_lcs.push(lc);
            if let Some(r) = quad_rsp {
                common_rsps.push(r);
            }
        } else {
            if matches!(
                config.bindata.output_type,
                BinDataOutputType::Spec | BinDataOutputType::Both
            ) {
                let pha_path = config
                    .output_dir
                    .join(format!("{prefix}_{gti_word}_clean_Q{q}.pha"));
                write_spectrum_file(&pha_path, &meta, &spec)?;
                summary.spectrum_paths.push(pha_path);

                if let Some(rsp) = quad_rsp {
                    let rsp_path = config.output_dir.join(&rsp_filename);
                    write_response_file(&rsp_path, &meta, &rsp)?;
                    summary.response_paths.push(rsp_path);
                }
            }

            if matches!(
                config.bindata.output_type,
                BinDataOutputType::Lc | BinDataOutputType::Both
            ) {
                let lc_path = config
                    .output_dir
                    .join(format!("{prefix}_{gti_word}_clean_Q{q}.lc"));
                write_lightcurve_file(&lc_path, &meta, &lc)?;
                summary.lightcurve_paths.push(lc_path);
            }
        }

        if let Some(we) = wevt {
            weighted_quads[q] = we;
            wrote_weighted = true;
        }
    }

    // Write combined all-quadrant `.pha`, `.lc`, `.rsp` in `GtiType::Common` mode
    if config.gti_type == GtiType::Common {
        let rsp_filename = if config.bindata.gen_resp && !common_rsps.is_empty() {
            format!("{prefix}_{gti_word}_clean.rsp")
        } else {
            "NONE".to_string()
        };
        if matches!(
            config.bindata.output_type,
            BinDataOutputType::Spec | BinDataOutputType::Both
        ) && !common_specs.is_empty()
        {
            let comb_spec = combine_spectra(&common_specs, false, rsp_filename.clone());
            let pha_path = config
                .output_dir
                .join(format!("{prefix}_{gti_word}_clean.pha"));
            write_spectrum_file(&pha_path, &meta, &comb_spec)?;
            summary.spectrum_paths.push(pha_path);

            if !common_rsps.is_empty() {
                let exps: Vec<f64> = common_specs.iter().map(|s| s.exposure).collect();
                let comb_rsp = combine_responses(&common_rsps, &exps, false);
                let rsp_path = config.output_dir.join(&rsp_filename);
                write_response_file(&rsp_path, &meta, &comb_rsp)?;
                summary.response_paths.push(rsp_path);
            }
        }
        if matches!(
            config.bindata.output_type,
            BinDataOutputType::Lc | BinDataOutputType::Both
        ) && !common_lcs.is_empty()
        {
            let comb_lc = combine_lightcurves(&common_lcs);
            let lc_path = config
                .output_dir
                .join(format!("{prefix}_{gti_word}_clean.lc"));
            write_lightcurve_file(&lc_path, &meta, &comb_lc)?;
            summary.lightcurve_paths.push(lc_path);
        }
    }

    if wrote_weighted {
        let wevt_path = config
            .output_dir
            .join(format!("{prefix}_{gti_word}_clean_weight.evt"));
        write_l2_event_file(
            &wevt_path,
            &meta,
            &weighted_quads,
            &combined_veto,
            &combined_ssm,
            &combined_temp,
            &final_gtis,
            &EventFileWriteOpts {
                gti_type: Some(config.gti_type),
                ncompat: Some(ncompat),
                gaincor_applied: true,
                lld_cut: config.evtencor.apply_lld,
                uld_cut: config.evtencor.apply_uld,
                pixel_exposure: Some(final_pix_exp.clone()),
            },
        )?;
        summary.weighted_evt_path = Some(wevt_path);
    }

    Ok(summary)
}

/// Standalone helper to run Stage 1 `science2event` directly from L1 FITS + TCT -> L2 FITS + bunch + hdr.
pub fn run_science2event(
    l1_file: impl AsRef<Path>,
    tct_file: impl AsRef<Path>,
    out_evt: impl AsRef<Path>,
    out_bunch: impl AsRef<Path>,
    out_hdr: impl AsRef<Path>,
    cfg: &Science2EventConfig,
) -> Result<[usize; NUM_QUADS]> {
    let l1_path = l1_file.as_ref();
    let meta = CztiObsMeta::from_fits(l1_path)?;
    let tct = CztiTct::from_fits(tct_file)?;
    let (gtis, _) = read_czti_gtis(l1_path)?;

    let decoded: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| -> Result<_> {
            let raw = read_l1_quad_packets(l1_path, q)?;
            Ok(decode_quad_packets(q as u8, &raw, cfg, |czt_sec| {
                tct.interpolate(czt_sec, CZTI_DEFAULT_FEB_LAG, None)
                    .unwrap_or(czt_sec)
            }))
        })
        .collect::<Result<Vec<_>>>()?;

    let counts = [
        decoded[0].events.len(),
        decoded[1].events.len(),
        decoded[2].events.len(),
        decoded[3].events.len(),
    ];

    let mut veto = VetoSpectrumTable::default();
    let mut ssm = SsmDataTable::default();
    let mut temp = TempTable::default();
    let mut hdr = FrameHeaderTable::default();
    for d in &decoded {
        veto.time.extend_from_slice(&d.veto_spec.time);
        veto.cztseccnt.extend_from_slice(&d.veto_spec.cztseccnt);
        veto.veto_spec.extend_from_slice(&d.veto_spec.veto_spec);
        veto.quad_id.extend_from_slice(&d.veto_spec.quad_id);
        ssm.time.extend_from_slice(&d.ssm_data.time);
        ssm.cztseccnt.extend_from_slice(&d.ssm_data.cztseccnt);
        ssm.veto_spec.extend_from_slice(&d.ssm_data.veto_spec);
        ssm.czt_spec.extend_from_slice(&d.ssm_data.czt_spec);
        ssm.czt_spec_veto.extend_from_slice(&d.ssm_data.czt_spec_veto);
        ssm.czt_spec_alpha.extend_from_slice(&d.ssm_data.czt_spec_alpha);
        ssm.quad_id.extend_from_slice(&d.ssm_data.quad_id);
        temp.time.extend_from_slice(&d.temp.time);
        temp.cztseccnt.extend_from_slice(&d.temp.cztseccnt);
        temp.temperature.extend_from_slice(&d.temp.temperature);
        temp.quad_flag.extend_from_slice(&d.temp.quad_flag);
        hdr.time.extend_from_slice(&d.headers.time);
        hdr.cztseccnt.extend_from_slice(&d.headers.cztseccnt);
        hdr.data_id.extend_from_slice(&d.headers.data_id);
        hdr.mode_id.extend_from_slice(&d.headers.mode_id);
        hdr.word_count.extend_from_slice(&d.headers.word_count);
        hdr.frame_no.extend_from_slice(&d.headers.frame_no);
        hdr.status_flag.extend_from_slice(&d.headers.status_flag);
        hdr.write_pkt_no.extend_from_slice(&d.headers.write_pkt_no);
        hdr.read_pkt_no.extend_from_slice(&d.headers.read_pkt_no);
        hdr.command.extend_from_slice(&d.headers.command);
        hdr.laxpc_time.extend_from_slice(&d.headers.laxpc_time);
        hdr.dcnt.extend_from_slice(&d.headers.dcnt);
        hdr.error_count.extend_from_slice(&d.headers.error_count);
        hdr.error_flag.extend_from_slice(&d.headers.error_flag);
        hdr.boot_page_no.extend_from_slice(&d.headers.boot_page_no);
        hdr.czt_no.extend_from_slice(&d.headers.czt_no);
        hdr.czt_status.extend_from_slice(&d.headers.czt_status);
        hdr.evt_read_mode.extend_from_slice(&d.headers.evt_read_mode);
        hdr.cmd_status.extend_from_slice(&d.headers.cmd_status);
        hdr.buffer_no.extend_from_slice(&d.headers.buffer_no);
        hdr.base_add.extend_from_slice(&d.headers.base_add);
        hdr.veto_spec_range.extend_from_slice(&d.headers.veto_spec_range);
        hdr.channel_no.extend_from_slice(&d.headers.channel_no);
        hdr.adc_output.extend_from_slice(&d.headers.adc_output);
        hdr.cmd_1sec.extend_from_slice(&d.headers.cmd_1sec);
        hdr.alpha_count.extend_from_slice(&d.headers.alpha_count);
        hdr.veto_count.extend_from_slice(&d.headers.veto_count);
        hdr.czt_count_lt_uld.extend_from_slice(&d.headers.czt_count_lt_uld);
        hdr.czt_count_gt_uld.extend_from_slice(&d.headers.czt_count_gt_uld);
        hdr.czt_data_read.extend_from_slice(&d.headers.czt_data_read);
    }

    let quads = [
        decoded[0].events.clone(),
        decoded[1].events.clone(),
        decoded[2].events.clone(),
        decoded[3].events.clone(),
    ];
    let bunches = [
        decoded[0].bunches.clone(),
        decoded[1].bunches.clone(),
        decoded[2].bunches.clone(),
        decoded[3].bunches.clone(),
    ];

    write_l2_event_file(
        out_evt,
        &meta,
        &quads,
        &veto,
        &ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts::default(),
    )?;
    write_bunch_file(out_bunch, &meta, &bunches)?;
    write_hdr_file(out_hdr, &meta, &hdr)?;

    Ok(counts)
}
