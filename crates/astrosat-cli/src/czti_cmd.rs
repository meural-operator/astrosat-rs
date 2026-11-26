//! `astrosat czti` — end-to-end and per-stage CZTI Level-1 -> Level-2 commands.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use astrosat_core::{CztiTct, MkfTable, CZTI_DEFAULT_FEB_LAG};
use astrosat_czti::bindata::{BinDataConfig, BinDataOutputType};
use astrosat_czti::bunchclean::{bunchclean_quad, BunchCleanConfig};
use astrosat_czti::datasel::datasel_quad;
use astrosat_czti::dpigen::{generate_dph_dpi, DpiGenConfig};
use astrosat_czti::evtclean::{evtclean_quad, EvtCleanConfig};
use astrosat_czti::evtencor::{evtencor_quad, EvtEnCorConfig};
use astrosat_czti::fits_io::{
    extract_mkf_rpy_rad, read_aux_tables, read_badpix_fits, read_czti_gtis, read_l1_quad_packets,
    read_l2_event_file, read_temp_extension, write_aspect_file, write_badpix_file,
    write_bunch_file, write_dph_file, write_dpi_file, write_gti_file, write_l2_event_file,
    write_livetime_file, write_sky_image_file, CztiCaldb, CztiObsMeta, EventFileWriteOpts,
};
use astrosat_czti::gaas::compute_quad_aspect;
use astrosat_czti::gtigen::{
    generate_czti_gtis, intervals_from_mkf_mask, parse_mkf_threshold_text, parse_user_gti_text,
    subtract_intervals, GtiGenConfig,
};
use astrosat_czti::image::{
    balanced_cross_correlation, balanced_cross_correlation_from_mask, combine_sky_images,
    fft_cross_correlation, ImageAlgorithm,
};
use astrosat_czti::noiseclean::{
    eventsep_quad, flickpixclean_quad_with_gti, heavybunchclean_quad_with_exposure,
    noisypixclean_quad, superbunchclean_quad, NoiseCleanConfig,
};
use astrosat_czti::pha2energy::{pha2energy_quad, QuadTempSeries};
use astrosat_czti::pipeline::{
    run_pipeline, run_science2event, CztiPipelineConfig, CztiPipelineSummary,
};
use astrosat_czti::pixclean::{pixclean_quad, PixCleanConfig};
use astrosat_czti::science2event::{decode_quad_packets, Science2EventConfig};
use astrosat_czti::stages::{run_addspec, run_bindata};
use astrosat_czti::types::{
    pix_flag, BadpixTable, BunchTable, EventSelectionAlgo, GtiType, QuadAspect, QuadEvents,
    QuadLivetime, TempTable, NUM_QUADS, PIXELS_PER_QUAD, QUAD_DIM,
};
use astrosat_fits::FitsFile;
use clap::ArgMatches;
use rayon::prelude::*;
use serde::Serialize;

use crate::cli::{
    is_cli_arg, resolve_flag, CztiAddspecArgs, CztiAlgoArg, CztiBindataArgs, CztiBunchcleanArgs,
    CztiCleanArgs, CztiCommand, CztiDataselArgs, CztiDecodeArgs, CztiDpigenArgs, CztiEvtcleanArgs,
    CztiEvtencorArgs, CztiGaasArgs, CztiGtiTypeArg, CztiGtigenArgs, CztiImageAlgoArg,
    CztiImageArgs, CztiOutputTypeArg, CztiPha2energyArgs, CztiRunArgs,
};
use crate::config::CztiFileConfig;
use crate::output::OutputMode;

#[derive(Debug, Serialize)]
struct CztiStageCountsReport {
    stage: &'static str,
    output: PathBuf,
    quad_events: [usize; NUM_QUADS],
    extra_outputs: Vec<PathBuf>,
}

#[derive(Debug, Serialize)]
struct CztiGtiReport {
    output: PathBuf,
    common_intervals: usize,
    common_exposure_sec: f64,
    quad_exposures_sec: [f64; NUM_QUADS],
}

#[derive(Debug, Serialize)]
struct CztiGaasReport {
    aspect_paths: Vec<PathBuf>,
    quad_aspect_deg: [(f64, f64, f64); NUM_QUADS],
}

#[derive(Debug, Serialize)]
struct CztiDpigenReport {
    out_dph: PathBuf,
    out_dpi: PathBuf,
    quad_total_counts: [i64; NUM_QUADS],
}

#[derive(Debug, Serialize)]
struct CztiImageReport {
    algorithm: String,
    image_paths: Vec<PathBuf>,
    peak_snr: Vec<f32>,
}

#[derive(Debug, Serialize)]
struct CztiBindataReport {
    ra_deg: f64,
    dec_deg: f64,
    outputs: Vec<PathBuf>,
}

#[derive(Debug, Serialize)]
struct CztiAddspecReport {
    inputs: usize,
    out_pha: PathBuf,
    out_rsp: Option<PathBuf>,
    exposure_sec: f64,
}

pub fn dispatch_czti(
    cmd: &CztiCommand,
    file_cfg: &CztiFileConfig,
    matches: &ArgMatches,
    out: OutputMode,
) -> Result<()> {
    match cmd {
        CztiCommand::Run(args) => {
            let res = run_czti_full(args, file_cfg, matches)?;
            out.emit(&res, print_czti_run_human)?;
        }
        CztiCommand::Decode(args) => {
            let res = run_czti_decode(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI Decode: Q0..Q3 events={:?} -> {}",
                    r.quad_events,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Bunchclean(args) => {
            let res = run_czti_bunchclean(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI Bunchclean: Q0..Q3 events={:?} -> {}",
                    r.quad_events,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Pha2energy(args) => {
            let res = run_czti_pha2energy(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI PHA2Energy: Q0..Q3 events={:?} -> {}",
                    r.quad_events,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Gtigen(args) => {
            let res = run_czti_gtigen(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI GTIGen: {} common intervals ({:.2} s), quad exp={:.2?} -> {}",
                    r.common_intervals,
                    r.common_exposure_sec,
                    r.quad_exposures_sec,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Gaas(args) => {
            let res = run_czti_gaas(args)?;
            out.emit(&res, |r| {
                println!("CZTI GAAS Aspect Computed:");
                for (q, (ra, dec, tw)) in r.quad_aspect_deg.iter().enumerate() {
                    println!("  Q{q}: RA={ra:.5} deg, DEC={dec:.5} deg, TWIST={tw:.5} deg");
                }
            })?;
        }
        CztiCommand::Datasel(args) => {
            let res = run_czti_datasel(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI Datasel: Q0..Q3 events={:?} -> {}",
                    r.quad_events,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Clean(args) => {
            let res = run_czti_clean(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI Clean: Q0..Q3 events={:?} -> {}",
                    r.quad_events,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Evtencor(args) => {
            let res = run_czti_evtencor(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI EvtEnCor: Q0..Q3 events={:?} -> {}",
                    r.quad_events,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Evtclean(args) => {
            let res = run_czti_evtclean(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI EvtClean: Q0..Q3 events={:?} -> {}",
                    r.quad_events,
                    r.output.display()
                );
            })?;
        }
        CztiCommand::Dpigen(args) => {
            let res = run_czti_dpigen(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI DPIGen: counts={:?} -> DPH={}, DPI={}",
                    r.quad_total_counts,
                    r.out_dph.display(),
                    r.out_dpi.display()
                );
            })?;
        }
        CztiCommand::Image(args) => {
            let res = run_czti_image(args)?;
            out.emit(&res, |r| {
                println!("CZTI Sky Image ({}):", r.algorithm);
                for (p, snr) in r.image_paths.iter().zip(r.peak_snr.iter()) {
                    println!("  {} (peak={:.2})", p.display(), snr);
                }
            })?;
        }
        CztiCommand::Bindata(args) => {
            let res = run_czti_bindata(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI Bindata (RA={:.4}, DEC={:.4}): {} output file(s)",
                    r.ra_deg,
                    r.dec_deg,
                    r.outputs.len(),
                );
            })?;
        }
        CztiCommand::Addspec(args) => {
            let res = run_czti_addspec(args)?;
            out.emit(&res, |r| {
                println!(
                    "CZTI Addspec ({} inputs, exposure={:.2} s) -> {}",
                    r.inputs,
                    r.exposure_sec,
                    r.out_pha.display()
                );
            })?;
        }
    }
    Ok(())
}

fn run_czti_full(
    args: &CztiRunArgs,
    fc: &CztiFileConfig,
    matches: &ArgMatches,
) -> Result<CztiPipelineSummary> {
    let input_dir = args.input.clone().or_else(|| fc.input_dir.clone());
    let discovered = input_dir.as_deref().map(discover_czti_l1).transpose()?;

    let m0 = args
        .m0
        .clone()
        .or_else(|| fc.l1_m0_file.clone())
        .or_else(|| discovered.as_ref().and_then(|d| d.m0.clone()))
        .context("missing CZTI modeM0 L1 FITS; pass --m0 <file> or --input <czti_l1_dir>")?;
    let ss = args
        .ss
        .clone()
        .or_else(|| fc.l1_ss_file.clone())
        .or_else(|| discovered.as_ref().and_then(|d| d.ss.clone()));
    let tct = args
        .tct
        .clone()
        .or_else(|| fc.tct_file.clone())
        .or_else(|| discovered.as_ref().and_then(|d| d.tct.clone()))
        .context("missing CZTI .tct file; pass --tct <file> or --input <czti_l1_dir>")?;
    let mkf = args
        .mkf
        .clone()
        .or_else(|| fc.mkf_file.clone())
        .or_else(|| discovered.as_ref().and_then(|d| d.mkf.clone()))
        .context("missing CZTI .mkf file; pass --mkf <file> or --input <czti_l1_dir>")?;
    let caldb = args
        .caldb
        .clone()
        .or_else(|| fc.caldb_dir.clone())
        .context("missing CZTI CALDB directory; pass --caldb <dir>")?;
    let output_dir = args
        .output
        .clone()
        .or_else(|| fc.output_dir.clone())
        .unwrap_or_else(|| PathBuf::from("czti_l2_out"));

    let mut cfg = CztiPipelineConfig::new(m0, ss, tct, mkf, caldb, output_dir);
    cfg.output_prefix = args.prefix.clone().or_else(|| fc.output_prefix.clone());
    cfg.mkf_threshold_file = args
        .mkf_thresholds
        .clone()
        .or_else(|| fc.mkf_threshold_file.clone())
        .or_else(|| discovered.as_ref().and_then(|d| d.mkf_thresholds.clone()));
    cfg.user_gti_file = args.user_gti.clone().or_else(|| fc.user_gti_file.clone());
    cfg.entry_stage = if is_cli_arg(matches, "entry_stage") {
        args.entry_stage
    } else {
        fc.entry_stage.unwrap_or(args.entry_stage)
    };
    cfg.exit_stage = if is_cli_arg(matches, "exit_stage") {
        args.exit_stage
    } else {
        fc.exit_stage.unwrap_or(args.exit_stage)
    };

    if let Some(gt) = args.gti_type {
        cfg.gti_type = match gt {
            CztiGtiTypeArg::Quad => GtiType::Quad,
            CztiGtiTypeArg::Common => GtiType::Common,
        };
    } else if let Some(ref s) = fc.gti_type {
        cfg.gti_type = GtiType::parse(s).context("invalid czti.gti_type in config")?;
    }

    if let Some(al) = args.algo {
        cfg.evt_sel_algo = match al {
            CztiAlgoArg::Default => EventSelectionAlgo::Default,
            CztiAlgoArg::Generalized => EventSelectionAlgo::Generalized,
        };
    } else if let Some(ref s) = fc.algo {
        cfg.evt_sel_algo = if s.eq_ignore_ascii_case("generalized") {
            EventSelectionAlgo::Generalized
        } else if s.eq_ignore_ascii_case("default") {
            EventSelectionAlgo::Default
        } else {
            bail!(
                "invalid czti.algo '{}' in config (expected 'default' or 'generalized')",
                s
            );
        };
    }

    let ra_opt = args.ra.or(fc.ra_deg);
    let dec_opt = args.dec.or(fc.dec_deg);
    match (ra_opt, dec_opt) {
        (Some(ra), Some(dec)) => {
            cfg.bindata.ra_src_deg = ra;
            cfg.bindata.dec_src_deg = dec;
            cfg.bindata.source_coords = Some((ra, dec));
        }
        (None, None) => {}
        _ => bail!("both --ra and --dec must be provided together (via CLI or config)"),
    }
    if let Some(tb) = args.tbin.or(fc.time_binsize) {
        cfg.bindata.time_binsize = tb;
    }
    if let Some(emin) = args.emin.or(fc.emin_kev) {
        cfg.bindata.emin_kev = emin;
    }
    if let Some(emax) = args.emax.or(fc.emax_kev) {
        cfg.bindata.emax_kev = emax;
    }
    cfg.bindata.output_type = if is_cli_arg(matches, "output_type") {
        match args.output_type {
            CztiOutputTypeArg::Spec => BinDataOutputType::Spec,
            CztiOutputTypeArg::Lc => BinDataOutputType::Lc,
            CztiOutputTypeArg::Both => BinDataOutputType::Both,
        }
    } else if let Some(ref ot) = fc.output_type {
        BinDataOutputType::parse(ot).context("invalid czti.output_type in config")?
    } else {
        BinDataOutputType::Both
    };
    cfg.bindata.mask_weight = resolve_flag(
        matches,
        "mask_weight",
        "no_mask_weight",
        fc.mask_weight,
        true,
    );
    cfg.bindata.gen_resp = resolve_flag(matches, "resp", "no_resp", fc.gen_resp, true);
    cfg.bindata.generate_eventfile = resolve_flag(
        matches,
        "gen_evt",
        "no_gen_evt",
        fc.generate_eventfile,
        false,
    );
    if let Some(ia) = args.image_algo {
        cfg.image.algorithm = match ia {
            CztiImageAlgoArg::Balcc => ImageAlgorithm::BalCC,
            CztiImageAlgoArg::Fft => ImageAlgorithm::Fft,
        };
    } else if let Some(ref s) = fc.image_algo {
        cfg.image.algorithm = ImageAlgorithm::parse(s).context("invalid czti.image_algo")?;
    }
    if let Some(ov) = args.oversampling.or(fc.oversampling) {
        if ov == 0 || ov > 64 {
            bail!("--oversampling must be in 1..=64, got {ov}");
        }
        cfg.image.oversampling_factor = ov;
    }
    if let Some(ref qs) = args.quads.clone().or_else(|| fc.quads.clone()) {
        cfg.image.quads = qs.clone();
        cfg.bindata.quads = qs.clone();
    }

    run_pipeline(&cfg).map_err(|e| anyhow::anyhow!("{e}"))
}

#[derive(Default)]
struct DiscoveredCzti {
    m0: Option<PathBuf>,
    ss: Option<PathBuf>,
    tct: Option<PathBuf>,
    mkf: Option<PathBuf>,
    mkf_thresholds: Option<PathBuf>,
}

fn discover_czti_l1(root: &Path) -> Result<DiscoveredCzti> {
    if !root.exists() {
        bail!("CZTI input directory does not exist: {}", root.display());
    }
    let mut files = Vec::new();
    collect_files_recursive(root, &mut files)?;
    files.sort();

    let mut d = DiscoveredCzti::default();
    for p in files {
        let name = p
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        if d.m0.is_none() && name.contains("cztm0") && name.ends_with(".fits") {
            d.m0 = Some(p);
        } else if d.ss.is_none() && name.contains("cztss") && name.ends_with(".fits") {
            d.ss = Some(p);
        } else if d.tct.is_none() && name.ends_with(".tct") {
            d.tct = Some(p);
        } else if d.mkf.is_none() && name.ends_with(".mkf") {
            d.mkf = Some(p);
        } else if d.mkf_thresholds.is_none() && name == "mkfthresholds.txt" {
            d.mkf_thresholds = Some(p);
        }
    }
    Ok(d)
}

fn collect_files_recursive(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let p = entry.path();
        if p.is_dir() {
            collect_files_recursive(&p, out)?;
        } else if p.is_file() {
            out.push(p);
        }
    }
    Ok(())
}

fn run_czti_decode(args: &CztiDecodeArgs) -> Result<CztiStageCountsReport> {
    let stem = args
        .output
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "czti".to_string());
    let dir = args
        .output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(dir)?;

    let bunch_path = args
        .bunch
        .clone()
        .unwrap_or_else(|| dir.join(format!("{stem}_bunch.fits")));
    let hdr_path = args
        .hdr
        .clone()
        .unwrap_or_else(|| dir.join(format!("{stem}.hdr")));

    let cfg = Science2EventConfig {
        big_endian: !args.little_endian,
    };
    let counts = run_science2event(
        &args.m0,
        &args.tct,
        &args.output,
        &bunch_path,
        &hdr_path,
        &cfg,
    )?;

    Ok(CztiStageCountsReport {
        stage: "decode",
        output: args.output.clone(),
        quad_events: counts,
        extra_outputs: vec![bunch_path, hdr_path],
    })
}

fn run_czti_bunchclean(args: &CztiBunchcleanArgs) -> Result<CztiStageCountsReport> {
    let (meta, quads, gtis, _pix_exp, gaincor_applied) = read_l2_event_file(&args.input)?;
    let (veto, ssm, temp) = read_aux_tables(&args.input);
    let bunches = read_bunch_file(&args.bunch)?;
    let cfg = if args.generalized {
        BunchCleanConfig::generalized_compat()
    } else {
        BunchCleanConfig::default()
    };

    let res: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| bunchclean_quad(&quads[q], &bunches[q], &gtis.quad[q], &cfg))
        .collect();

    let ncompat = res[0].ncompat;
    let out_quads: [QuadEvents; NUM_QUADS] = [
        res[0].events.clone(),
        res[1].events.clone(),
        res[2].events.clone(),
        res[3].events.clone(),
    ];
    let out_bunches: [BunchTable; NUM_QUADS] = [
        res[0].out_bunches.clone(),
        res[1].out_bunches.clone(),
        res[2].out_bunches.clone(),
        res[3].out_bunches.clone(),
    ];
    let out_lv: [QuadLivetime; NUM_QUADS] = [
        res[0].livetime.clone(),
        res[1].livetime.clone(),
        res[2].livetime.clone(),
        res[3].livetime.clone(),
    ];
    let pix_exp: [Vec<f64>; NUM_QUADS] = [
        res[0].pixel_exposure.clone(),
        res[1].pixel_exposure.clone(),
        res[2].pixel_exposure.clone(),
        res[3].pixel_exposure.clone(),
    ];
    let counts = [
        out_quads[0].len(),
        out_quads[1].len(),
        out_quads[2].len(),
        out_quads[3].len(),
    ];

    ensure_parent_dir(&args.output)?;
    write_l2_event_file(
        &args.output,
        &meta,
        &out_quads,
        &veto,
        &ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts {
            ncompat: Some(ncompat),
            gaincor_applied,
            pixel_exposure: Some(pix_exp),
            ..Default::default()
        },
    )?;

    let mut extra = Vec::new();
    if let Some(ref bp) = args.out_bunch {
        ensure_parent_dir(bp)?;
        write_bunch_file(bp, &meta, &out_bunches)?;
        extra.push(bp.clone());
    }
    if let Some(ref lp) = args.out_livetime {
        ensure_parent_dir(lp)?;
        write_livetime_file(lp, &meta, &out_lv)?;
        extra.push(lp.clone());
    }

    Ok(CztiStageCountsReport {
        stage: "bunchclean",
        output: args.output.clone(),
        quad_events: counts,
        extra_outputs: extra,
    })
}

fn run_czti_pha2energy(args: &CztiPha2energyArgs) -> Result<CztiStageCountsReport> {
    let (meta, mut quads, gtis, pix_exp, gaincor_applied) = read_l2_event_file(&args.input)?;
    let (veto, ssm, mut temp) = read_aux_tables(&args.input);
    if let Some(ref ss_path) = args.ss {
        if !ss_path.is_file() {
            bail!("--ss file does not exist: {}", ss_path.display());
        }
        if let Ok(t) = read_temp_extension(ss_path) {
            temp = t;
        } else {
            // Check if `--ss` is a Level-1 modeSS packet file (`CZT_QUAD1..4`).
            let ss_fits = FitsFile::open(ss_path)
                .with_context(|| format!("failed to open --ss FITS {}", ss_path.display()))?;
            if ss_fits.hdu_by_name("CZT_QUAD1").is_err() {
                bail!(
                    "--ss {} contains neither a TEMP extension nor Level-1 CZT_QUAD1 telemetry",
                    ss_path.display()
                );
            }
            let tct_path = args.tct.as_ref().with_context(|| {
                format!(
                    "--ss {} is a Level-1 modeSS packet file; pass --tct <tct_file> to decode TEMP",
                    ss_path.display()
                )
            })?;
            let tct = CztiTct::from_fits(tct_path)?;
            let cfg_s2e = Science2EventConfig::default();
            let decoded_ss: Vec<_> = (0..NUM_QUADS)
                .into_par_iter()
                .map(|q| {
                    let raw = read_l1_quad_packets(ss_path, q)?;
                    Ok(decode_quad_packets(q as u8, &raw, &cfg_s2e, |czt_sec| {
                        tct.interpolate(czt_sec, CZTI_DEFAULT_FEB_LAG, None)
                            .unwrap_or(czt_sec)
                    }))
                })
                .collect::<astrosat_czti::Result<Vec<_>>>()?;
            let mut decoded_temp = TempTable::default();
            for s in decoded_ss {
                decoded_temp.time.extend(s.temp.time);
                decoded_temp.cztseccnt.extend(s.temp.cztseccnt);
                decoded_temp.temperature.extend(s.temp.temperature);
                decoded_temp.quad_flag.extend(s.temp.quad_flag);
            }
            if decoded_temp.time.is_empty() {
                bail!(
                    "decoded 0 TEMP rows from Level-1 modeSS file {}",
                    ss_path.display()
                );
            }
            temp = decoded_temp;
        }
    }

    let caldb = CztiCaldb::open(&args.caldb)?;
    let ebounds = caldb.read_ebounds(meta.tstart)?;
    let gain_offsets: Vec<_> = (0..NUM_QUADS)
        .map(|q| caldb.read_gain_offset(q, meta.tstart))
        .collect::<Result<Vec<_>, _>>()?;
    let temps: Vec<QuadTempSeries> = (0..NUM_QUADS)
        .map(|q| QuadTempSeries::from_temp_table(&temp, q as u8))
        .collect();

    quads
        .par_iter_mut()
        .enumerate()
        .for_each(|(q, ev)| pha2energy_quad(ev, &temps[q], &gain_offsets[q], &ebounds));

    let counts = [quads[0].len(), quads[1].len(), quads[2].len(), quads[3].len()];
    ensure_parent_dir(&args.output)?;
    write_l2_event_file(
        &args.output,
        &meta,
        &quads,
        &veto,
        &ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts {
            gaincor_applied,
            pixel_exposure: pix_exp,
            ..Default::default()
        },
    )?;

    Ok(CztiStageCountsReport {
        stage: "pha2energy",
        output: args.output.clone(),
        quad_events: counts,
        extra_outputs: Vec::new(),
    })
}

fn run_czti_gtigen(args: &CztiGtigenArgs) -> Result<CztiGtiReport> {
    let meta = CztiObsMeta::from_fits(&args.m0)?;
    let (m0_gtis, m0_bti) = read_czti_gtis(&args.m0)?;
    let mkf = MkfTable::from_fits(&args.mkf)?;
    let rules = if let Some(ref tf) = args.mkf_thresholds {
        let txt = fs::read_to_string(tf)?;
        parse_mkf_threshold_text(&txt)
    } else {
        GtiGenConfig::default().thresholds
    };

    let mut mask = vec![true; mkf.len()];
    for rule in &rules {
        if !rule.enabled {
            continue;
        }
        if let Some(col) = mkf.column(&rule.column) {
            for (i, &val) in col.iter().enumerate() {
                if !val.is_finite()
                    || rule.min.is_some_and(|lo| val < lo)
                    || rule.max.is_some_and(|hi| val > hi)
                {
                    mask[i] = false;
                }
            }
        }
    }
    let mkf_gti = intervals_from_mkf_mask(mkf.times(), &mask);
    let user_gti = args
        .user_gti
        .as_deref()
        .map(|p| -> Result<_> {
            let txt = fs::read_to_string(p)?;
            Ok(parse_user_gti_text(&txt)?)
        })
        .transpose()?;

    let final_gtis = generate_czti_gtis(
        &mkf_gti,
        &m0_gtis.common,
        &m0_gtis.quad,
        if m0_bti.is_empty() { None } else { Some(&m0_bti) },
        user_gti.as_deref(),
        args.min_exposure,
    )?;

    ensure_parent_dir(&args.output)?;
    write_gti_file(&args.output, &meta, &final_gtis)?;

    let sum_exp = |ivs: &[(f64, f64)]| -> f64 { ivs.iter().map(|&(a, b)| (b - a).max(0.0)).sum() };

    Ok(CztiGtiReport {
        output: args.output.clone(),
        common_intervals: final_gtis.common.len(),
        common_exposure_sec: sum_exp(&final_gtis.common),
        quad_exposures_sec: [
            sum_exp(&final_gtis.quad[0]),
            sum_exp(&final_gtis.quad[1]),
            sum_exp(&final_gtis.quad[2]),
            sum_exp(&final_gtis.quad[3]),
        ],
    })
}

fn run_czti_gaas(args: &CztiGaasArgs) -> Result<CztiGaasReport> {
    let mkf = MkfTable::from_fits(&args.mkf)?;
    let rpy = extract_mkf_rpy_rad(&mkf)?;
    let t0 = mkf.times().first().copied().unwrap_or(0.0);
    let t1 = mkf.times().last().copied().unwrap_or(0.0);
    let caldb = CztiCaldb::open(&args.caldb)?;
    let gtis = args
        .gti
        .as_deref()
        .map(read_czti_gtis)
        .transpose()?
        .map(|(g, _)| g);

    let mut meta = CztiObsMeta::from_fits(&args.mkf).unwrap_or_default();
    meta.tstart = t0;
    meta.tstop = t1;

    fs::create_dir_all(&args.output_dir)?;
    let mut aspect_paths = Vec::with_capacity(NUM_QUADS);
    let mut quad_aspect_deg = [(0.0, 0.0, 0.0); NUM_QUADS];

    for q in 0..NUM_QUADS {
        let (tname, tmat) = caldb.read_teldef(q, t0)?;
        let qgti = gtis.as_ref().and_then(|g| match args.gti_type {
            CztiGtiTypeArg::Common => (!g.common.is_empty()).then_some(g.common.as_slice()),
            CztiGtiTypeArg::Quad => (!g.quad[q].is_empty()).then_some(g.quad[q].as_slice()),
        });
        let asp = compute_quad_aspect(q as u8, tname, &tmat, &rpy, t0, t1, qgti);
        quad_aspect_deg[q] = (asp.ra_deg, asp.dec_deg, asp.twist_deg);
        let p = args.output_dir.join(format!("{}.aspect_Q{q}", args.prefix));
        write_aspect_file(&p, &meta, &asp)?;
        aspect_paths.push(p);
    }

    Ok(CztiGaasReport {
        aspect_paths,
        quad_aspect_deg,
    })
}

fn run_czti_datasel(args: &CztiDataselArgs) -> Result<CztiStageCountsReport> {
    let (meta, quads, _old_gtis, pix_exp, gaincor_applied) = read_l2_event_file(&args.input)?;
    let (veto, ssm, temp) = read_aux_tables(&args.input);
    let (gtis, _) = read_czti_gtis(&args.gti)?;

    let out_quads: Vec<QuadEvents> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| {
            let qgti = match args.gti_type {
                CztiGtiTypeArg::Common => &gtis.common,
                CztiGtiTypeArg::Quad => &gtis.quad[q],
            };
            datasel_quad(&quads[q], qgti)
        })
        .collect();
    let arr = into_quad_array(out_quads);
    let counts = [arr[0].len(), arr[1].len(), arr[2].len(), arr[3].len()];

    ensure_parent_dir(&args.output)?;
    write_l2_event_file(
        &args.output,
        &meta,
        &arr,
        &veto,
        &ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts {
            gti_type: Some(match args.gti_type {
                CztiGtiTypeArg::Quad => GtiType::Quad,
                CztiGtiTypeArg::Common => GtiType::Common,
            }),
            gaincor_applied,
            pixel_exposure: pix_exp,
            ..Default::default()
        },
    )?;

    Ok(CztiStageCountsReport {
        stage: "datasel",
        output: args.output.clone(),
        quad_events: counts,
        extra_outputs: Vec::new(),
    })
}

fn run_czti_clean(args: &CztiCleanArgs) -> Result<CztiStageCountsReport> {
    let (meta, quads, mut gtis, pix_exp_opt, gaincor_applied) = read_l2_event_file(&args.input)?;
    let (veto, ssm, temp) = read_aux_tables(&args.input);
    let caldb = CztiCaldb::open(&args.caldb)?;
    let caldb_bp = caldb.read_badpix(meta.tstart)?;
    let pix_exp = pix_exp_opt.unwrap_or_else(|| {
        [
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
        ]
    });
    let livetimes = args
        .livetime
        .as_deref()
        .map(read_livetime_file)
        .transpose()?
        .unwrap_or_else(|| default_livetimes(&quads));
    let bunches = args
        .bunch
        .as_deref()
        .map(read_bunch_file)
        .transpose()?
        .unwrap_or_default();

    let pc_cfg = PixCleanConfig::default();
    let nc_cfg = NoiseCleanConfig::default();

    let per_q: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| {
            let qgti = match args.gti_type {
                CztiGtiTypeArg::Common => &gtis.common,
                CztiGtiTypeArg::Quad => &gtis.quad[q],
            };
            match args.algo {
                CztiAlgoArg::Default => {
                    let pc = pixclean_quad(
                        &quads[q],
                        &livetimes[q],
                        qgti,
                        &caldb_bp.flags[q],
                        &pix_exp[q],
                        &pc_cfg,
                    );
                    (
                        pc.single_events,
                        pc.livetime,
                        pc.badpix_flags,
                        pc.pixel_exposure,
                        Vec::new(),
                    )
                }
                CztiAlgoArg::Generalized => {
                    let gc = caldb.read_gaincor_lld_uld(q, meta.tstart).ok();
                    let lld = gc.map(|g| g.lld).unwrap_or_else(|| vec![0i16; PIXELS_PER_QUAD]);
                    let np = noisypixclean_quad(&quads[q], &caldb_bp.flags[q], &pix_exp[q], &nc_cfg);
                    let sb = superbunchclean_quad(&np.events, &bunches[q], &livetimes[q], &nc_cfg);
                    let (hb_ev, hb_pexp) = heavybunchclean_quad_with_exposure(
                        &sb.events,
                        &bunches[q],
                        &lld,
                        &np.pixel_exposure,
                        &nc_cfg,
                    );
                    let fp = flickpixclean_quad_with_gti(
                        &hb_ev,
                        &np.badpix_flags,
                        &hb_pexp,
                        Some(qgti),
                        &nc_cfg,
                    );
                    let (mut singles, _) = eventsep_quad(&fp.events, nc_cfg.double_time_sec);
                    if !sb.bti_intervals.is_empty() {
                        let updated_gti = subtract_intervals(qgti, &sb.bti_intervals);
                        singles.exposure =
                            updated_gti.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
                    }
                    (
                        singles,
                        sb.livetime,
                        fp.badpix_flags,
                        fp.pixel_exposure,
                        sb.bti_intervals,
                    )
                }
            }
        })
        .collect();

    let mut out_q = Vec::with_capacity(NUM_QUADS);
    let mut out_lv = Vec::with_capacity(NUM_QUADS);
    let mut out_bp = BadpixTable::default();
    let mut out_pe: [Vec<f64>; NUM_QUADS] = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    for (q, (ev, lv, bp, pe, bti)) in per_q.into_iter().enumerate() {
        out_q.push(ev);
        out_lv.push(lv);
        out_bp.flags[q] = bp;
        out_pe[q] = pe;
        if !bti.is_empty() {
            gtis.quad[q] = subtract_intervals(&gtis.quad[q], &bti);
        }
    }

    let q_arr = into_quad_array(out_q);
    let lv_arr = into_quad_array(out_lv);
    let counts = [q_arr[0].len(), q_arr[1].len(), q_arr[2].len(), q_arr[3].len()];

    ensure_parent_dir(&args.output)?;
    write_l2_event_file(
        &args.output,
        &meta,
        &q_arr,
        &veto,
        &ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts {
            ncompat: matches!(args.algo, CztiAlgoArg::Generalized).then_some(1),
            gaincor_applied,
            pixel_exposure: Some(out_pe),
            ..Default::default()
        },
    )?;

    let mut extra = Vec::new();
    if let Some(ref bp_path) = args.out_badpix {
        ensure_parent_dir(bp_path)?;
        write_badpix_file(bp_path, &meta, &out_bp)?;
        extra.push(bp_path.clone());
    }
    if let Some(ref lv_path) = args.out_livetime {
        ensure_parent_dir(lv_path)?;
        write_livetime_file(lv_path, &meta, &lv_arr)?;
        extra.push(lv_path.clone());
    }

    Ok(CztiStageCountsReport {
        stage: "clean",
        output: args.output.clone(),
        quad_events: counts,
        extra_outputs: extra,
    })
}

fn run_czti_evtencor(args: &CztiEvtencorArgs) -> Result<CztiStageCountsReport> {
    let (meta, quads, gtis, pix_exp, gaincor_applied) = read_l2_event_file(&args.input)?;
    let (veto, ssm, temp) = read_aux_tables(&args.input);
    let caldb = CztiCaldb::open(&args.caldb)?;
    let gc: Vec<_> = (0..NUM_QUADS)
        .map(|q| caldb.read_gaincor_lld_uld(q, meta.tstart))
        .collect::<Result<Vec<_>, _>>()?;

    let cfg = EvtEnCorConfig {
        apply_lld: !args.no_lld,
        apply_uld: !args.no_uld,
        time_gcor: !args.no_time_gcor,
        rng_seed: match args.seed {
            Some(s) => {
                let b = s.clamp(1, i32::MAX - 12);
                (b, b + 6, b + 12)
            }
            None => (5, 11, 17),
        },
    };

    let out_q: Vec<QuadEvents> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| evtencor_quad(&quads[q], &gc[q], meta.tstart, gaincor_applied, &cfg).events)
        .collect();
    let arr = into_quad_array(out_q);
    let counts = [arr[0].len(), arr[1].len(), arr[2].len(), arr[3].len()];

    ensure_parent_dir(&args.output)?;
    write_l2_event_file(
        &args.output,
        &meta,
        &arr,
        &veto,
        &ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts {
            gaincor_applied: true,
            lld_cut: cfg.apply_lld,
            uld_cut: cfg.apply_uld,
            pixel_exposure: pix_exp,
            ..Default::default()
        },
    )?;

    Ok(CztiStageCountsReport {
        stage: "evtencor",
        output: args.output.clone(),
        quad_events: counts,
        extra_outputs: Vec::new(),
    })
}

fn run_czti_evtclean(args: &CztiEvtcleanArgs) -> Result<CztiStageCountsReport> {
    let (meta, quads, gtis, pix_exp, gaincor_applied) = read_l2_event_file(&args.input)?;
    let (veto, ssm, temp) = read_aux_tables(&args.input);
    let badpix = args.badpix.as_deref().map(read_badpix_fits).transpose()?;

    let cfg = EvtCleanConfig {
        allowed_alpha: vec![args.alpha],
        veto_range: (args.veto_min, args.veto_max),
        is_double_event: false,
    };

    let out_q: Vec<QuadEvents> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| {
            let ev = evtclean_quad(&quads[q], &cfg);
            if let Some(ref bp) = badpix {
                let keep: Vec<bool> = (0..ev.len())
                    .map(|i| bp.get(q, ev.detid[i], ev.pixid[i]) < pix_flag::NOISY)
                    .collect();
                ev.filter_by_mask(&keep)
            } else {
                ev
            }
        })
        .collect();
    let arr = into_quad_array(out_q);
    let counts = [arr[0].len(), arr[1].len(), arr[2].len(), arr[3].len()];

    ensure_parent_dir(&args.output)?;
    write_l2_event_file(
        &args.output,
        &meta,
        &arr,
        &veto,
        &ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts {
            gaincor_applied,
            pixel_exposure: pix_exp,
            ..Default::default()
        },
    )?;

    Ok(CztiStageCountsReport {
        stage: "evtclean",
        output: args.output.clone(),
        quad_events: counts,
        extra_outputs: Vec::new(),
    })
}

fn run_czti_dpigen(args: &CztiDpigenArgs) -> Result<CztiDpigenReport> {
    let (meta, quads, _gtis, _pix_exp, _) = read_l2_event_file(&args.input)?;
    let caldb = CztiCaldb::open(&args.caldb)?;
    let badpix = read_badpix_fits(&args.badpix)?;
    let eff = [
        caldb.read_eff_area(0, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(1, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(2, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(3, meta.tstart).unwrap_or_default(),
    ];
    let cfg = DpiGenConfig {
        energy_range: Some((args.emin, args.emax)),
        ..Default::default()
    };
    let prod = generate_dph_dpi(&quads, &badpix, &eff, &cfg);
    let gt = match args.gti_type {
        CztiGtiTypeArg::Quad => GtiType::Quad,
        CztiGtiTypeArg::Common => GtiType::Common,
    };

    ensure_parent_dir(&args.out_dph)?;
    ensure_parent_dir(&args.out_dpi)?;
    write_dph_file(&args.out_dph, &meta, gt, &prod.quad_dph, &prod.full_dph)?;
    write_dpi_file(&args.out_dpi, &meta, gt, &prod.quad_dpi, &prod.full_dpi)?;

    let sum_q = |q: usize| -> i64 {
        prod.quad_dph[q]
            .iter()
            .flat_map(|r| r.iter())
            .map(|&v| v as i64)
            .sum()
    };

    Ok(CztiDpigenReport {
        out_dph: args.out_dph.clone(),
        out_dpi: args.out_dpi.clone(),
        quad_total_counts: [sum_q(0), sum_q(1), sum_q(2), sum_q(3)],
    })
}

fn run_czti_image(args: &CztiImageArgs) -> Result<CztiImageReport> {
    if args.oversampling == 0 || args.oversampling > 64 {
        bail!("--oversampling must be in 1..=64, got {}", args.oversampling);
    }
    let meta = CztiObsMeta::from_fits(&args.dpi)?;
    let dpi_fits = FitsFile::open(&args.dpi)?;
    let caldb = CztiCaldb::open(&args.caldb)?;
    fs::create_dir_all(&args.output_dir)?;

    let mut image_paths = Vec::new();
    let mut peak_snr = Vec::new();
    let algo_name = match args.algo {
        CztiImageAlgoArg::Balcc => "balcc",
        CztiImageAlgoArg::Fft => "fft",
    };

    let shadowlib_all = if matches!(args.algo, CztiImageAlgoArg::Balcc) {
        caldb.read_shadowlib_all(meta.tstart).ok()
    } else {
        None
    };

    let mut sky_imgs = Vec::new();
    let mut first_aspect = None;
    for &q in &args.quads {
        if q >= NUM_QUADS {
            continue;
        }
        let hdu = dpi_fits
            .hdu_by_name(&format!("Q{q}"))
            .or_else(|_| dpi_fits.hdu(q + 1))?;
        let img = hdu.as_image()?;
        let flat = img.read_f32()?;
        let mut dpi_2d = vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM];
        for y in 0..QUAD_DIM {
            for x in 0..QUAD_DIM {
                if y * QUAD_DIM + x < flat.len() {
                    dpi_2d[y][x] = flat[y * QUAD_DIM + x];
                }
            }
        }
        let asp = read_quad_aspect(&args.aspect, q)?;
        if first_aspect.is_none() {
            first_aspect = Some(asp.clone());
        }
        let sky_img = match args.algo {
            CztiImageAlgoArg::Balcc => {
                if let Some((grid_dim, ref all_shadows)) = shadowlib_all {
                    let refs: Vec<&[f32]> = all_shadows[q].iter().map(|v| v.as_slice()).collect();
                    balanced_cross_correlation(&dpi_2d, &refs, grid_dim, &asp)
                } else {
                    let mask = caldb
                        .read_mask_pattern(q, meta.tstart)
                        .unwrap_or_else(|_| vec![vec![0.5; QUAD_DIM]; QUAD_DIM]);
                    balanced_cross_correlation_from_mask(&dpi_2d, &mask, 15, &asp)
                }
            }
            CztiImageAlgoArg::Fft => {
                let mask = caldb
                    .read_mask_pattern(q, meta.tstart)
                    .unwrap_or_else(|_| vec![vec![0.5; QUAD_DIM]; QUAD_DIM]);
                fft_cross_correlation(&dpi_2d, &mask, args.oversampling, &asp)
            }
        };

        let pk = sky_img
            .data
            .iter()
            .copied()
            .fold(f32::NEG_INFINITY, f32::max);
        peak_snr.push(pk);
        let p = args.output_dir.join(format!("{}.img_Q{q}", args.prefix));
        write_sky_image_file(&p, &format!("Q{q}"), &meta, &sky_img)?;
        image_paths.push(p);
        sky_imgs.push(sky_img);
    }

    if sky_imgs.len() > 1 {
        if let Some(ref asp0) = first_aspect {
            let combined = combine_sky_images(&sky_imgs, asp0);
            let pk = combined
                .data
                .iter()
                .copied()
                .fold(f32::NEG_INFINITY, f32::max);
            peak_snr.push(pk);
            let p = args.output_dir.join(format!("{}.img_all", args.prefix));
            write_sky_image_file(&p, "IMAGE", &meta, &combined)?;
            image_paths.push(p);
        }
    }

    Ok(CztiImageReport {
        algorithm: algo_name.to_string(),
        image_paths,
        peak_snr,
    })
}

fn run_czti_bindata(args: &CztiBindataArgs) -> Result<CztiBindataReport> {
    fs::create_dir_all(&args.output_dir)?;
    let out_prefix = args.output_dir.join(&args.prefix);
    let source_coords = match (args.ra, args.dec) {
        (Some(r), Some(d)) => Some((r, d)),
        (None, None) => None,
        _ => bail!("both --ra and --dec must be provided together"),
    };
    let cfg = BinDataConfig {
        ra_src_deg: args.ra.unwrap_or(0.0),
        dec_src_deg: args.dec.unwrap_or(0.0),
        source_coords,
        mask_weight: !args.no_mask_weight,
        output_type: match args.output_type {
            CztiOutputTypeArg::Spec => BinDataOutputType::Spec,
            CztiOutputTypeArg::Lc => BinDataOutputType::Lc,
            CztiOutputTypeArg::Both => BinDataOutputType::Both,
        },
        time_binsize: args.tbin,
        emin_kev: args.emin,
        emax_kev: args.emax,
        gen_resp: !args.no_resp,
        generate_eventfile: args.gen_evt,
        gti_type: match args.gti_type {
            CztiGtiTypeArg::Quad => GtiType::Quad,
            CztiGtiTypeArg::Common => GtiType::Common,
        },
        quads: args.quads.clone(),
        ..Default::default()
    };

    let outputs = run_bindata(
        &args.input,
        &args.mkf,
        &args.badpix,
        &args.livetime,
        &args.caldb,
        &out_prefix,
        &cfg,
    )?;

    let (ra_pnt, dec_pnt) = read_mkf_pointing_radec(&args.mkf);
    let (ra_deg, dec_deg) = cfg.resolve_source_radec(ra_pnt, dec_pnt);

    Ok(CztiBindataReport {
        ra_deg,
        dec_deg,
        outputs,
    })
}

fn read_mkf_pointing_radec(mkf_path: &Path) -> (f64, f64) {
    let Ok(fits) = FitsFile::open(mkf_path) else {
        return (0.0, 0.0);
    };
    let hdu = fits
        .hdu_by_name("MKF")
        .or_else(|_| fits.hdu_by_name("MKF_DATA"))
        .or_else(|_| fits.hdu(1))
        .or_else(|_| fits.hdu(0));
    let ra_pnt = hdu
        .as_ref()
        .ok()
        .and_then(|h| h.header().get_f64("RA_PNT"))
        .or_else(|| fits.hdu(0).ok().and_then(|h| h.header().get_f64("RA_PNT")))
        .unwrap_or(0.0);
    let dec_pnt = hdu
        .as_ref()
        .ok()
        .and_then(|h| h.header().get_f64("DEC_PNT"))
        .or_else(|| fits.hdu(0).ok().and_then(|h| h.header().get_f64("DEC_PNT")))
        .unwrap_or(0.0);
    (ra_pnt, dec_pnt)
}

fn run_czti_addspec(args: &CztiAddspecArgs) -> Result<CztiAddspecReport> {
    ensure_parent_dir(&args.out_pha)?;
    if let Some(ref rsp_out) = args.out_rsp {
        ensure_parent_dir(rsp_out)?;
    }
    let (spec, _rsp) = run_addspec(
        &args.pha,
        &args.rsp,
        &args.out_pha,
        args.out_rsp.as_deref(),
        args.add_exposure,
    )?;
    Ok(CztiAddspecReport {
        inputs: args.pha.len(),
        out_pha: args.out_pha.clone(),
        out_rsp: args.out_rsp.clone(),
        exposure_sec: spec.exposure,
    })
}

fn read_bunch_file(path: &Path) -> Result<[BunchTable; NUM_QUADS]> {
    let fits = FitsFile::open(path)?;
    let mut vec = Vec::with_capacity(NUM_QUADS);
    for q in 0..NUM_QUADS {
        let h = fits.hdu_by_name(&format!("Q{q}")).or_else(|_| fits.hdu(q + 1))?;
        let bt = h.as_bintable()?;
        vec.push(BunchTable {
            quad_id: q as u8,
            time: bt.read_col_f64("TIME")?,
            time_dfs: bt.read_col_u8("TIME_DFS")?,
            time_dsl: bt.read_col_u8("TIME_DSL")?,
            num_event: bt.read_col_u8("NUMEVENT")?,
            detid1: bt.read_col_u8("DETID1")?,
            detid2: bt.read_col_u8("DETID2")?,
            detid3: bt.read_col_u8("DETID3")?,
            detid4: bt.read_col_u8("DETID4")?,
            row_num: bt.read_col_i32("ROWNUM")?,
            detid_fevt: bt.read_col_u8("DETID_FEVT")?,
            pixid_fevt: bt.read_col_u8("PIXID_FEVT")?,
            detid_sevt: bt.read_col_u8("DETID_SEVT")?,
            pixid_sevt: bt.read_col_u8("PIXID_SEVT")?,
            detid_levt: bt.read_col_u8("DETID_LEVT")?,
            pixid_levt: bt.read_col_u8("PIXID_LEVT")?,
        });
    }
    Ok(into_quad_array(vec))
}

fn read_livetime_file(path: &Path) -> Result<[QuadLivetime; NUM_QUADS]> {
    let fits = FitsFile::open(path)?;
    let mut vec = Vec::with_capacity(NUM_QUADS);
    for q in 0..NUM_QUADS {
        let h = fits.hdu_by_name(&format!("Q{q}")).or_else(|_| fits.hdu(q + 1))?;
        let bt = h.as_bintable()?;
        let binsize = h.header().get_f64("LV_BINSIZE").unwrap_or(1.0);
        vec.push(QuadLivetime {
            quad_id: q as u8,
            time: bt.read_col_f64("TIME")?,
            fracexp: bt.read_col_f64("FRACEXP")?,
            binsize,
        });
    }
    Ok(into_quad_array(vec))
}

fn default_livetimes(quads: &[QuadEvents; NUM_QUADS]) -> [QuadLivetime; NUM_QUADS] {
    let mk = |q: usize| {
        let t0 = quads[q].tstart.floor();
        let t1 = quads[q].tstop.ceil().max(t0 + 1.0);
        let n = ((t1 - t0) as usize).max(1);
        QuadLivetime {
            quad_id: q as u8,
            time: (0..n).map(|i| t0 + (i as f64) + 0.5).collect(),
            fracexp: vec![1.0; n],
            binsize: 1.0,
        }
    };
    [mk(0), mk(1), mk(2), mk(3)]
}

fn read_quad_aspect(path_or_dir: &Path, q: usize) -> Result<QuadAspect> {
    let target = if path_or_dir.is_dir() {
        let mut found = None;
        for entry in fs::read_dir(path_or_dir)? {
            let p = entry?.path();
            let s = p.to_string_lossy();
            if s.ends_with(&format!("aspect_Q{q}")) {
                found = Some(p);
                break;
            }
        }
        found.with_context(|| {
            format!("no *aspect_Q{q} file found in {}", path_or_dir.display())
        })?
    } else if q > 0 {
        let s = path_or_dir.to_string_lossy();
        if s.contains("aspect_Q0") {
            PathBuf::from(s.replace("aspect_Q0", &format!("aspect_Q{q}")))
        } else {
            path_or_dir.to_path_buf()
        }
    } else {
        path_or_dir.to_path_buf()
    };

    let fits = FitsFile::open(&target)?;
    let h = fits.hdu_by_name("ASPECT").or_else(|_| fits.hdu(1))?;
    let hdr = h.header();
    let bt = h.as_bintable()?;
    Ok(QuadAspect {
        quad_id: q as u8,
        teldef_name: hdr.get_str("TELDEF").unwrap_or("TELDEF").trim().to_string(),
        time: bt.read_col_f64("TIME").unwrap_or_default(),
        nx: bt.read_col_f32("NX").unwrap_or_default(),
        ny: bt.read_col_f32("NY").unwrap_or_default(),
        nz: bt.read_col_f32("NZ").unwrap_or_default(),
        nxt: bt.read_col_f32("NXT").unwrap_or_default(),
        nyt: bt.read_col_f32("NYT").unwrap_or_default(),
        nzt: bt.read_col_f32("NZT").unwrap_or_default(),
        ra_deg: hdr.get_f64("RA").unwrap_or(0.0),
        dec_deg: hdr.get_f64("DEC").unwrap_or(0.0),
        twist_deg: hdr.get_f64("TWIST").unwrap_or(0.0),
    })
}

fn into_quad_array<T>(mut v: Vec<T>) -> [T; NUM_QUADS] {
    assert_eq!(v.len(), NUM_QUADS);
    let mut it = v.drain(..);
    [
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
    ]
}

fn ensure_parent_dir(path: &Path) -> Result<()> {
    if let Some(p) = path.parent() {
        if !p.as_os_str().is_empty() {
            fs::create_dir_all(p)?;
        }
    }
    Ok(())
}

fn print_czti_run_human(s: &CztiPipelineSummary) {
    println!(
        "CZTI Level-1 -> Level-2 Pipeline Complete (OBS_ID: {})",
        s.obs_id.as_deref().unwrap_or("-")
    );
    println!("  Stage 1 Raw Events        : {:?}", s.stage1_raw_events);
    println!("  Stage 1 Bunchclean Events : {:?}", s.stage1_bunchclean_events);
    println!("  Stage 2 Clean Events      : {:?}", s.stage2_clean_events);
    println!("  Quad Exposures (s)        : {:.2?}", s.quad_exposures_sec);
    if let Some(ref p) = s.clean_evt_path {
        println!("  Clean Event FITS          : {}", p.display());
    }
    if let Some(ref p) = s.gti_path {
        println!("  GTI FITS                  : {}", p.display());
    }
    if let Some(ref p) = s.dpi_path {
        println!("  DPI FITS                  : {}", p.display());
    }
    println!(
        "  Products                  : {} images, {} spectra, {} light curves, {} responses",
        s.image_paths.len(),
        s.spectrum_paths.len(),
        s.lightcurve_paths.len(),
        s.response_paths.len()
    );
}
