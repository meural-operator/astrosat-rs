//! `astrosat laxpc` — end-to-end and per-stage LAXPC Level-1 -> Level-2 commands.

use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{bail, Context, Result};
use astrosat_core::{GtiInterval, GtiSet};
use astrosat_fits::{FitsFile, OgipMeta};
use astrosat_laxpc::background::{
    background_lightcurve_shift_fraction, compute_gain_shift, shift_background_spectrum,
    BackgroundCatalog, BackgroundModel2D, BackgroundPha, ChannelMap, GainShiftSolution,
    GainTable, HvJumpTable,
};
use astrosat_laxpc::bb::BbcSeries;
use astrosat_laxpc::ea::decode_ea_orbits;
use astrosat_laxpc::frames::{FrameDiagnostics, FrameMode};
use astrosat_laxpc::gti::{
    generate_standard_gti, load_laxpc_tct_files, load_merged_mkf_files, GtiConfig,
};
use astrosat_laxpc::pipeline::{
    read_and_sort_l1_fits_frames, run_pipeline, write_orbit_from_mkf, MultiUnitPipelineResult,
    PipelineConfig,
};
use astrosat_laxpc::products::{
    bin_events_raw, build_calibrated_events, build_lightcurve, build_spectrum,
    compute_gti_frame_gap_loss, write_ascii_lc, write_ascii_pha, write_format_b_event_fits,
    write_gti_fits, write_lightcurve_fits, write_spectrum_fits,
};
use astrosat_laxpc::types::{utc_to_mission_day, AnodeSelection, EventFilter, LaxpcUnit};
use clap::ArgMatches;
use serde::Serialize;

use crate::cli::{
    resolve_flag, LaxpcBackshiftArgs, LaxpcCommand, LaxpcDecodeArgs, LaxpcGtigenArgs,
    LaxpcLightcurveArgs, LaxpcOrbitArgs, LaxpcRunArgs, LaxpcSpectrumArgs,
};
use crate::config::LaxpcFileConfig;
use crate::output::OutputMode;

#[derive(Debug, Serialize)]
struct LaxpcDecodeReport {
    unit: String,
    ea_files: usize,
    bb_files: usize,
    ea_diagnostics: FrameDiagnostics,
    bb_diagnostics: Option<FrameDiagnostics>,
    total_valid_events: usize,
    frame_gaps: usize,
    output_evn: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct LaxpcGtigenReport {
    unit: String,
    output: PathBuf,
    intervals: usize,
    gross_exposure_sec: f64,
    frame_gap_loss_sec: f64,
    net_exposure_sec: f64,
    tct_rejected_glitches: usize,
}

#[derive(Debug, Serialize)]
struct LaxpcLightcurveReport {
    unit: String,
    output: PathBuf,
    ascii_output: Option<PathBuf>,
    tbin: f64,
    bins: usize,
    mean_rate_cps: f64,
    exposure_sec: f64,
}

#[derive(Debug, Serialize)]
struct LaxpcSpectrumReport {
    unit: String,
    output: PathBuf,
    ascii_output: Option<PathBuf>,
    channels: usize,
    exposure_sec: f64,
    total_corrected_rate_cps: f64,
}

#[derive(Debug, Serialize)]
struct LaxpcBackshiftReport {
    unit: String,
    mission_day: f64,
    gain_shift: GainShiftSolution,
    lcbk_shift_fraction: Option<f64>,
    shifted_spectrum_fits: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct LaxpcOrbitReport {
    mkf_files: usize,
    output: PathBuf,
    tstart: f64,
    tstop: f64,
}

pub fn dispatch_laxpc(
    cmd: &LaxpcCommand,
    file_cfg: &LaxpcFileConfig,
    matches: &ArgMatches,
    out: OutputMode,
) -> Result<()> {
    match cmd {
        LaxpcCommand::Run(args) => {
            let res = run_laxpc_full(args, file_cfg, matches)?;
            out.emit(&res, print_laxpc_run_human)?;
        }
        LaxpcCommand::Decode(args) => {
            let res = run_laxpc_decode(args)?;
            out.emit(&res, |r| {
                println!(
                    "LAXPC Decode [{}]: {} valid events across {} retained EA frames ({} gaps)",
                    r.unit,
                    r.total_valid_events,
                    r.ea_diagnostics.retained_frames,
                    r.frame_gaps
                );
                if let Some(ref p) = r.output_evn {
                    println!("  Wrote Format B events: {}", p.display());
                }
            })?;
        }
        LaxpcCommand::Gtigen(args) => {
            let res = run_laxpc_gtigen(args)?;
            out.emit(&res, |r| {
                println!(
                    "LAXPC GTI [{}]: {} intervals, net exposure={:.3} s (gross={:.3} s, gap_loss={:.3} s) -> {}",
                    r.unit,
                    r.intervals,
                    r.net_exposure_sec,
                    r.gross_exposure_sec,
                    r.frame_gap_loss_sec,
                    r.output.display()
                );
            })?;
        }
        LaxpcCommand::Lightcurve(args) => {
            let res = run_laxpc_lightcurve(args)?;
            out.emit(&res, |r| {
                println!(
                    "LAXPC Light Curve [{}]: {} bins (tbin={} s), mean={:.3} cps, exp={:.3} s -> {}",
                    r.unit,
                    r.bins,
                    r.tbin,
                    r.mean_rate_cps,
                    r.exposure_sec,
                    r.output.display()
                );
            })?;
        }
        LaxpcCommand::Spectrum(args) => {
            let res = run_laxpc_spectrum(args)?;
            out.emit(&res, |r| {
                println!(
                    "LAXPC Spectrum [{}]: {} channels, rate={:.3} cps, exp={:.3} s -> {}",
                    r.unit,
                    r.channels,
                    r.total_corrected_rate_cps,
                    r.exposure_sec,
                    r.output.display()
                );
            })?;
        }
        LaxpcCommand::Backshift(args) => {
            let res = run_laxpc_backshift(args)?;
            out.emit(&res, |r| {
                println!(
                    "LAXPC Backshift [{}] (mission day {:.2}): delta_ch_30keV={:.3}, fsh={:.5}, RMF={}",
                    r.unit,
                    r.mission_day,
                    r.gain_shift.delta_ch_30kev,
                    r.gain_shift.secular_scale,
                    r.gain_shift.recommended_rmf
                );
                if let Some(fsh) = r.lcbk_shift_fraction {
                    println!("  LC background shift fshbk: {fsh:.6}");
                }
                if let Some(ref p) = r.shifted_spectrum_fits {
                    println!("  Wrote shifted background spectrum: {}", p.display());
                }
            })?;
        }
        LaxpcCommand::Orbit(args) => {
            let res = run_laxpc_orbit(args)?;
            out.emit(&res, |r| {
                println!(
                    "LAXPC Orbit: {} MKF file(s), T={:.2}..{:.2} -> {}",
                    r.mkf_files,
                    r.tstart,
                    r.tstop,
                    r.output.display()
                );
            })?;
        }
    }
    Ok(())
}

fn run_laxpc_full(
    args: &LaxpcRunArgs,
    fc: &LaxpcFileConfig,
    matches: &ArgMatches,
) -> Result<MultiUnitPipelineResult> {
    let input_dir = args.input.clone().or_else(|| fc.input_dir.clone());
    if input_dir.is_none() {
        bail!("missing LAXPC Level-1 input directory; pass --input <l1_dir> or set laxpc.input_dir in --config");
    }

    let units = if let Some(ref units_raw) = args.units.clone().or_else(|| fc.units.clone()) {
        parse_units(units_raw)?
    } else {
        LaxpcUnit::ALL.to_vec()
    };

    let caldb_dir = args.caldb.clone().or_else(|| fc.caldb_dir.clone());

    let layer1_only = resolve_flag(
        matches,
        "layer1_only",
        "no_layer1_only",
        fc.layer1_only,
        false,
    );
    let anodes = if let Some(ref a) = args.anodes {
        parse_anodes(a)?
    } else if crate::cli::is_cli_arg(matches, "layer1_only") {
        AnodeSelection::Layer1
    } else if crate::cli::is_cli_arg(matches, "no_layer1_only") {
        AnodeSelection::AllMain
    } else if let Some(ref a) = fc.anodes {
        parse_anodes(a)?
    } else if layer1_only {
        AnodeSelection::Layer1
    } else {
        AnodeSelection::AllMain
    };

    let energy_range_kev = if let Some(er) = args.energy_range.as_deref().or(fc.energy_range.as_deref()) {
        let range = parse_energy_range(er)?;
        validate_channel_maps(caldb_dir.as_deref(), &units)?;
        Some(range)
    } else {
        None
    };

    let pha_min = args.pha_min.or(fc.pha_min).unwrap_or(0);
    let pha_max = args.pha_max.or(fc.pha_max).unwrap_or(1023);
    validate_pha_range(pha_min, pha_max)?;

    let filter = EventFilter {
        single_only: resolve_flag(
            matches,
            "single_only",
            "no_single_only",
            fc.single_only,
            false,
        ),
        anodes,
        pha_min,
        pha_max,
        energy_range_kev,
    };

    let gti_defaults = GtiConfig::default();
    let gti_config = GtiConfig {
        min_elv_deg: args.elv_min.or(fc.min_elv_deg).unwrap_or(gti_defaults.min_elv_deg),
        saa_entry_lon_deg: args
            .saa_entry_lon
            .or(fc.saa_entry_lon_deg)
            .unwrap_or(gti_defaults.saa_entry_lon_deg),
        saa_exit_lon_deg: args
            .saa_exit_lon
            .or(fc.saa_exit_lon_deg)
            .unwrap_or(gti_defaults.saa_exit_lon_deg),
        merge_gap_sec: args
            .merge_gap
            .or(fc.merge_gap_sec)
            .unwrap_or(gti_defaults.merge_gap_sec),
        min_duration_sec: args
            .min_duration
            .or(fc.min_duration_sec)
            .unwrap_or(gti_defaults.min_duration_sec),
        ..gti_defaults
    };

    let user_gti_utc = args
        .user_gti
        .clone()
        .or_else(|| fc.user_gti_file.clone())
        .as_deref()
        .map(load_user_gti)
        .transpose()?;

    let defaults = PipelineConfig::default();
    let cfg = PipelineConfig {
        input_dir,
        output_dir: args
            .output
            .clone()
            .or_else(|| fc.output_dir.clone())
            .unwrap_or(defaults.output_dir),
        units,
        caldb_dir,
        tbin: args.tbin.or(fc.tbin).unwrap_or(defaults.tbin),
        deadtime_sec: args.deadtime.or(fc.deadtime_sec).unwrap_or(defaults.deadtime_sec),
        spectrum_channels: args.channels.or(fc.spectrum_channels),
        write_format_a_combined: resolve_flag(
            matches,
            "format_a",
            "no_format_a",
            fc.write_format_a_combined,
            false,
        ),
        write_format_b_events: resolve_flag(
            matches,
            "format_b",
            "no_format_b",
            fc.write_format_b_events,
            true,
        ),
        write_orbit_file: resolve_flag(
            matches,
            "orbit",
            "no_orbit",
            fc.write_orbit_file,
            true,
        ),
        write_ascii_products: resolve_flag(
            matches,
            "ascii",
            "no_ascii",
            fc.write_ascii_products,
            false,
        ),
        filter,
        gti_config,
        user_gti_utc,
        use_user_gti_only: resolve_flag(
            matches,
            "user_gti_only",
            "no_user_gti_only",
            fc.user_gti_only,
            false,
        ),
        background_model_path: args
            .bkg_model
            .clone()
            .or_else(|| fc.background_model_path.clone()),
        background_pha_path: args
            .bkg_pha
            .clone()
            .or_else(|| fc.background_pha_path.clone()),
        background_peak_ch: args.bkg_peak_ch.or(fc.background_peak_ch),
        ..defaults
    };

    run_pipeline(&cfg).map_err(|e| anyhow::anyhow!("{e}"))
}

fn run_laxpc_decode(args: &LaxpcDecodeArgs) -> Result<LaxpcDecodeReport> {
    let unit = resolve_laxpc_unit(args.unit.as_deref(), &args.ea)?;
    let mut ea_diag = FrameDiagnostics::default();
    let mut ea_batches = Vec::with_capacity(args.ea.len());
    for p in &args.ea {
        let (frames, d) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Ea)?;
        merge_diag(&mut ea_diag, d);
        ea_batches.push(frames);
    }

    let mut bb_diag_opt = None;
    if !args.bb.is_empty() {
        let mut bd = FrameDiagnostics::default();
        for p in &args.bb {
            let (_frames, d) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Bb)?;
            merge_diag(&mut bd, d);
        }
        bb_diag_opt = Some(bd);
    }

    let stream = decode_ea_orbits(&ea_batches, unit);

    let output_evn = if let Some(ref out_path) = args.output {
        if args.tct.is_empty() {
            bail!("writing Format B .evn output requires --tct <tct_files...>");
        }
        let tct = load_laxpc_tct_files(&args.tct)?;
        let gti_utc = if let Some(ref gp) = args.gti {
            load_user_gti(gp)?
        } else {
            let t0 = tct.interpolate(stream.t_ins_min);
            let t1 = tct.interpolate(stream.t_ins_max);
            GtiSet::from_pairs(&[(t0, t1)])
        };
        let gti_ins_pairs: Vec<(f64, f64)> = gti_utc
            .intervals()
            .iter()
            .map(|iv| {
                (
                    astrosat_laxpc::gti::utc_to_ins(&tct, iv.start),
                    astrosat_laxpc::gti::utc_to_ins(&tct, iv.stop),
                )
            })
            .collect();
        let gti_utc_pairs: Vec<(f64, f64)> = gti_utc
            .intervals()
            .iter()
            .map(|iv| (iv.start, iv.stop))
            .collect();

        let filter = EventFilter::default();
        let nch = unit.default_spectrum_channels();
        let cal = build_calibrated_events(&stream.events, &gti_ins_pairs, &filter, nch, |ti| {
            tct.interpolate(ti)
        });
        let mut meta = OgipMeta {
            instrume: "LAXPC".into(),
            detnam: Some(unit.detnam().into()),
            exposure: gti_utc.exposure(),
            tstart: gti_utc.tstart().unwrap_or(0.0),
            tstop: gti_utc.tstop().unwrap_or(0.0),
            ..Default::default()
        };
        if let Some(first_ea) = args.ea.first() {
            if let Ok(fits) = FitsFile::open(first_ea) {
                for idx in [0usize, 1usize] {
                    if let Ok(h) = fits.hdu(idx) {
                        let hdr = h.header();
                        if meta.obs_id.is_none() {
                            meta.obs_id = hdr.get_str("OBS_ID").map(|s| s.trim().into());
                        }
                        if meta.object.is_none() {
                            meta.object = hdr.get_str("OBJECT").map(|s| s.trim().into());
                        }
                        if meta.ra_obj.is_none() {
                            meta.ra_obj = hdr
                                .get_f64("RA_NOM")
                                .or_else(|| hdr.get_f64("RA_OBJ"))
                                .or_else(|| hdr.get_f64("RA_PNT"));
                        }
                        if meta.dec_obj.is_none() {
                            meta.dec_obj = hdr
                                .get_f64("DEC_NOM")
                                .or_else(|| hdr.get_f64("DEC_OBJ"))
                                .or_else(|| hdr.get_f64("DEC_PNT"));
                        }
                    }
                }
            }
        }
        ensure_parent_dir(out_path)?;
        write_format_b_event_fits(out_path, &cal, &gti_utc_pairs, unit, &filter, nch, &meta)?;
        Some(out_path.clone())
    } else {
        None
    };

    Ok(LaxpcDecodeReport {
        unit: unit.detnam().to_string(),
        ea_files: args.ea.len(),
        bb_files: args.bb.len(),
        ea_diagnostics: ea_diag,
        bb_diagnostics: bb_diag_opt,
        total_valid_events: stream.total_valid_events,
        frame_gaps: stream.frame_gaps.len(),
        output_evn,
    })
}

fn run_laxpc_gtigen(args: &LaxpcGtigenArgs) -> Result<LaxpcGtigenReport> {
    let unit = resolve_laxpc_unit(args.unit.as_deref(), &args.ea)?;
    let mut ea_batches = Vec::with_capacity(args.ea.len());
    for p in &args.ea {
        let (frames, _) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Ea)?;
        ea_batches.push(frames);
    }
    let bbc = if !args.bb.is_empty() {
        let mut bb_batches = Vec::with_capacity(args.bb.len());
        for p in &args.bb {
            let (frames, _) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Bb)?;
            bb_batches.push(frames);
        }
        Some(BbcSeries::from_ordered_orbits(&bb_batches))
    } else {
        None
    };

    let tct = load_laxpc_tct_files(&args.tct)?;
    let mkf = load_merged_mkf_files(&args.mkf)?;
    let stream = decode_ea_orbits(&ea_batches, unit);
    let filter = EventFilter::default();
    let raw_bins = bin_events_raw(&stream.events, 1.0, &filter)?;
    let user_gti = args.user_gti.as_deref().map(load_user_gti).transpose()?;

    let gti_cfg = GtiConfig {
        min_elv_deg: args.elv_min,
        ..Default::default()
    };

    let res = generate_standard_gti(
        unit,
        &raw_bins,
        &tct,
        &mkf,
        bbc.as_ref(),
        &stream.frame_gaps,
        &gti_cfg,
        user_gti.as_ref(),
    )?;

    let gti_pairs: Vec<(f64, f64)> = res
        .gti_utc
        .intervals()
        .iter()
        .map(|iv| (iv.start, iv.stop))
        .collect();
    let mut meta = read_mkf_meta(&args.mkf, unit);
    meta.tstart = res.gti_utc.tstart().unwrap_or(0.0);
    meta.tstop = res.gti_utc.tstop().unwrap_or(0.0);
    meta.exposure = res.net_exposure_sec;
    meta.ontime = Some(res.gross_exposure_sec);
    meta.livetime = Some(res.net_exposure_sec);

    ensure_parent_dir(&args.output)?;
    write_gti_fits(&args.output, &gti_pairs, &meta)?;

    Ok(LaxpcGtigenReport {
        unit: unit.detnam().to_string(),
        output: args.output.clone(),
        intervals: gti_pairs.len(),
        gross_exposure_sec: res.gross_exposure_sec,
        frame_gap_loss_sec: res.frame_gap_loss_sec,
        net_exposure_sec: res.net_exposure_sec,
        tct_rejected_glitches: res.tct_rejected_glitches,
    })
}

fn run_laxpc_lightcurve(args: &LaxpcLightcurveArgs) -> Result<LaxpcLightcurveReport> {
    let unit = resolve_laxpc_unit(args.unit.as_deref(), &args.ea)?;
    let mut ea_batches = Vec::with_capacity(args.ea.len());
    for p in &args.ea {
        let (frames, _) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Ea)?;
        ea_batches.push(frames);
    }
    let bbc = if !args.bb.is_empty() {
        let mut bb_batches = Vec::with_capacity(args.bb.len());
        for p in &args.bb {
            let (frames, _) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Bb)?;
            bb_batches.push(frames);
        }
        Some(BbcSeries::from_ordered_orbits(&bb_batches))
    } else {
        None
    };

    let tct = load_laxpc_tct_files(&args.tct)?;
    let mkf = load_merged_mkf_files(&args.mkf)?;
    let stream = decode_ea_orbits(&ea_batches, unit);

    let mut filter = EventFilter {
        single_only: args.single_only,
        anodes: if args.layer1_only {
            AnodeSelection::Layer1
        } else {
            AnodeSelection::AllMain
        },
        pha_min: args.pha_min.unwrap_or(0),
        pha_max: args.pha_max.unwrap_or(1023),
        ..Default::default()
    };
    validate_pha_range(filter.pha_min, filter.pha_max)?;

    if let Some(ref er) = args.energy_range {
        let (elo, ehi) = parse_energy_range(er)?;
        validate_channel_maps(args.caldb.as_deref(), &[unit])?;
        let cmap_path = args
            .caldb
            .as_ref()
            .unwrap()
            .join(format!("{}0.channel", unit.prefix()));
        let cmap = ChannelMap::from_file(&cmap_path)?;
        let (p_lo, p_hi) = cmap.energy_to_pha_1024(elo, ehi);
        filter.pha_min = p_lo;
        filter.pha_max = p_hi;
        filter.energy_range_kev = Some((elo, ehi));
    } else if let Some(ref caldb_dir) = args.caldb {
        if !caldb_dir.exists() {
            bail!("--caldb directory does not exist: {}", caldb_dir.display());
        }
    }

    let gti_tbin = if args.tbin < 0.2 { 1.0 } else { args.tbin };
    let gti_raw = bin_events_raw(&stream.events, gti_tbin, &filter)?;
    let user_gti = args.gti.as_deref().map(load_user_gti).transpose()?;
    let (gti_ins, gti_utc, net_exposure_sec) = if let Some(ugti) = user_gti {
        let gti_ins = GtiSet::new(
            ugti.intervals()
                .iter()
                .map(|iv| {
                    GtiInterval::new(
                        astrosat_laxpc::gti::utc_to_ins(&tct, iv.start),
                        astrosat_laxpc::gti::utc_to_ins(&tct, iv.stop),
                    )
                })
                .collect(),
        );
        let ins_pairs: Vec<(f64, f64)> = gti_ins.intervals().iter().map(|iv| (iv.start, iv.stop)).collect();
        let gap_loss = compute_gti_frame_gap_loss(&ins_pairs, &stream.frame_gaps);
        let net_exp = (gti_ins.exposure() - gap_loss).max(0.0);
        (gti_ins, ugti, net_exp)
    } else {
        let gcfg = GtiConfig {
            deadtime_sec: args.deadtime,
            ..Default::default()
        };
        let gr = generate_standard_gti(
            unit,
            &gti_raw,
            &tct,
            &mkf,
            bbc.as_ref(),
            &stream.frame_gaps,
            &gcfg,
            None,
        )?;
        (gr.gti_ins, gr.gti_utc, gr.net_exposure_sec)
    };

    let gti_ins_pairs: Vec<(f64, f64)> = gti_ins.intervals().iter().map(|iv| (iv.start, iv.stop)).collect();
    let gti_utc_pairs: Vec<(f64, f64)> = gti_utc.intervals().iter().map(|iv| (iv.start, iv.stop)).collect();

    let lc_raw = if (args.tbin - gti_tbin).abs() < 1e-12 {
        gti_raw
    } else {
        bin_events_raw(&stream.events, args.tbin, &filter)?
    };

    let bkg_model = args
        .bkg_model
        .as_deref()
        .map(BackgroundModel2D::from_file)
        .transpose()?;
    let bkg_per_bin: Option<Vec<f64>> = bkg_model.as_ref().map(|m| {
        let mkf_t = mkf.times();
        let lon_col = mkf.column_any(&["EARTHLON", "SAT_LON", "LONGITUDE", "LON"]);
        let lat_col = mkf.column_any(&["EARTHLAT", "SAT_LAT", "LATITUDE", "LAT"]);
        (0..lc_raw.all_main_counts.len())
            .map(|i| {
                let tc = lc_raw.t0_ins + (i as f64 + 0.5) * lc_raw.tbin;
                let xb = tct.interpolate(tc);
                let lon = lon_col
                    .map(|c| interp_lon_deg(mkf_t, c, xb))
                    .unwrap_or(100.0);
                let lat = lat_col
                    .map(|c| interp_linear(mkf_t, c, xb))
                    .unwrap_or(0.0);
                m.evaluate(lon, lat, xb)
            })
            .collect()
    });

    let lc = build_lightcurve(
        &lc_raw,
        &gti_ins_pairs,
        &stream.frame_gaps,
        unit,
        args.deadtime,
        |ti| tct.interpolate(ti),
        bkg_per_bin.as_deref(),
    );

    let mut meta = read_mkf_meta(&args.mkf, unit);
    meta.tstart = gti_utc.tstart().unwrap_or(0.0);
    meta.tstop = gti_utc.tstop().unwrap_or(0.0);
    meta.exposure = net_exposure_sec;

    ensure_parent_dir(&args.output)?;
    write_lightcurve_fits(
        &args.output,
        &lc,
        &lc.rate_cps,
        &lc.error_cps,
        &gti_utc_pairs,
        &filter,
        args.deadtime,
        false,
        &meta,
    )?;
    if let Some(ref asc) = args.ascii {
        ensure_parent_dir(asc)?;
        write_ascii_lc(asc, &lc)?;
    }

    Ok(LaxpcLightcurveReport {
        unit: unit.detnam().to_string(),
        output: args.output.clone(),
        ascii_output: args.ascii.clone(),
        tbin: lc.tbin,
        bins: lc.time_utc.len(),
        mean_rate_cps: lc.mean_rate_cps,
        exposure_sec: meta.exposure,
    })
}

fn run_laxpc_spectrum(args: &LaxpcSpectrumArgs) -> Result<LaxpcSpectrumReport> {
    let unit = resolve_laxpc_unit(args.unit.as_deref(), &args.ea)?;
    let mut ea_batches = Vec::with_capacity(args.ea.len());
    for p in &args.ea {
        let (frames, _) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Ea)?;
        ea_batches.push(frames);
    }
    let bbc = if !args.bb.is_empty() {
        let mut bb_batches = Vec::with_capacity(args.bb.len());
        for p in &args.bb {
            let (frames, _) = read_and_sort_l1_fits_frames(p, unit, FrameMode::Bb)?;
            bb_batches.push(frames);
        }
        Some(BbcSeries::from_ordered_orbits(&bb_batches))
    } else {
        None
    };

    let tct = load_laxpc_tct_files(&args.tct)?;
    let mkf = load_merged_mkf_files(&args.mkf)?;
    let stream = decode_ea_orbits(&ea_batches, unit);

    let mut filter = EventFilter {
        single_only: args.single_only,
        anodes: if args.layer1_only {
            AnodeSelection::Layer1
        } else {
            AnodeSelection::AllMain
        },
        pha_min: args.pha_min.unwrap_or(0),
        pha_max: args.pha_max.unwrap_or(1023),
        ..Default::default()
    };
    validate_pha_range(filter.pha_min, filter.pha_max)?;

    if let Some(ref er) = args.energy_range {
        let (elo, ehi) = parse_energy_range(er)?;
        validate_channel_maps(args.caldb.as_deref(), &[unit])?;
        let cmap_path = args
            .caldb
            .as_ref()
            .unwrap()
            .join(format!("{}0.channel", unit.prefix()));
        let cmap = ChannelMap::from_file(&cmap_path)?;
        let (p_lo, p_hi) = cmap.energy_to_pha_1024(elo, ehi);
        filter.pha_min = p_lo;
        filter.pha_max = p_hi;
        filter.energy_range_kev = Some((elo, ehi));
    } else if let Some(ref caldb_dir) = args.caldb {
        if !caldb_dir.exists() {
            bail!("--caldb directory does not exist: {}", caldb_dir.display());
        }
    }

    let user_gti = args.gti.as_deref().map(load_user_gti).transpose()?;
    let (gti_ins, gti_utc) = if let Some(ugti) = user_gti {
        let ins = GtiSet::new(
            ugti.intervals()
                .iter()
                .map(|iv| {
                    GtiInterval::new(
                        astrosat_laxpc::gti::utc_to_ins(&tct, iv.start),
                        astrosat_laxpc::gti::utc_to_ins(&tct, iv.stop),
                    )
                })
                .collect(),
        );
        (ins, ugti)
    } else {
        let gti_raw = bin_events_raw(&stream.events, 1.0, &filter)?;
        let gcfg = GtiConfig {
            deadtime_sec: args.deadtime,
            ..Default::default()
        };
        let gr = generate_standard_gti(
            unit,
            &gti_raw,
            &tct,
            &mkf,
            bbc.as_ref(),
            &stream.frame_gaps,
            &gcfg,
            None,
        )?;
        (gr.gti_ins, gr.gti_utc)
    };

    let gti_ins_pairs: Vec<(f64, f64)> = gti_ins.intervals().iter().map(|iv| (iv.start, iv.stop)).collect();
    let nch = args.channels.unwrap_or_else(|| unit.default_spectrum_channels());
    let mut spec = build_spectrum(
        &stream.events,
        &gti_ins_pairs,
        &stream.frame_gaps,
        unit,
        &filter,
        nch,
        args.deadtime,
    )?;

    if filter.pha_min > 0 || filter.pha_max < 1023 {
        let ncf = (1024 / nch).max(1);
        let ch_lo = (filter.pha_min as usize) / ncf;
        let ch_hi = (filter.pha_max as usize) / ncf;
        for ch in 0..nch {
            if ch < ch_lo || ch > ch_hi {
                spec.selected_rate[ch] = 0.0;
                spec.selected_stat_err[ch] = 0.0;
                spec.anode_rates[ch] = [0.0; 11];
            }
        }
        spec.total_corrected_rate_cps = spec.selected_rate.iter().sum();
    }

    let mut meta = read_mkf_meta(&args.mkf, unit);
    meta.tstart = gti_utc.tstart().unwrap_or(0.0);
    meta.tstop = gti_utc.tstop().unwrap_or(0.0);
    meta.exposure = spec.exposure_sec;

    ensure_parent_dir(&args.output)?;
    write_spectrum_fits(
        &args.output,
        unit,
        &spec.selected_rate,
        &spec.selected_stat_err,
        spec.exposure_sec,
        args.deadtime,
        &filter,
        "NONE",
        "NONE",
        &meta,
    )?;

    if let Some(ref asc) = args.ascii {
        ensure_parent_dir(asc)?;
        let mday = utc_to_mission_day(0.5 * (meta.tstart + meta.tstop)).round() as i32;
        write_ascii_pha(asc, &spec, mday, filter.single_only)?;
    }

    Ok(LaxpcSpectrumReport {
        unit: unit.detnam().to_string(),
        output: args.output.clone(),
        ascii_output: args.ascii.clone(),
        channels: nch,
        exposure_sec: spec.exposure_sec,
        total_corrected_rate_cps: spec.total_corrected_rate_cps,
    })
}

fn run_laxpc_backshift(args: &LaxpcBackshiftArgs) -> Result<LaxpcBackshiftReport> {
    let unit = LaxpcUnit::from_str(&args.unit).map_err(|e| anyhow::anyhow!("{e}"))?;
    let gain_path = args.caldb.join("gaina8.dat");
    let gain_table = GainTable::from_file(&gain_path)
        .with_context(|| format!("failed to load {}", gain_path.display()))?;

    let hv_table = [args.caldb.join("lxhv.dat"), args.caldb.join("lxhv.data")]
        .into_iter()
        .find(|p| p.is_file())
        .map(|p| HvJumpTable::from_file(&p))
        .transpose()?;

    let mission_day = if let Some(md) = args.mission_day {
        md
    } else if let Some(ref fp) = args.fits {
        let fits = FitsFile::open(fp)?;
        let ts = fits
            .hdu(1)
            .ok()
            .and_then(|h| h.header().get_f64("TSTART"))
            .or_else(|| fits.hdu(0).ok().and_then(|h| h.header().get_f64("TSTART")))
            .context("could not find TSTART in FITS file; pass --mission-day")?;
        utc_to_mission_day(ts).max(0.0)
    } else {
        bail!("provide either --mission-day <days> or --fits <l2_fits>");
    };

    let bkg_pha = args
        .bkg_pha
        .as_deref()
        .map(BackgroundPha::from_file)
        .transpose()?;
    let bkg_day = bkg_pha.as_ref().map(|b| b.mission_day).unwrap_or(mission_day);
    let bkg_peak_ch = args.bkg_peak_ch.or_else(|| {
        if bkg_pha.is_some() {
            let cat = BackgroundCatalog::load_or_builtin(&[
                Some(args.caldb.as_path()),
                args.bkg_pha.as_deref().and_then(Path::parent),
            ]);
            cat.lookup_peak_30kev(unit, args.bkg_pha.as_deref(), bkg_day)
        } else {
            None
        }
    });

    let gain_shift = compute_gain_shift(
        unit,
        mission_day,
        bkg_day,
        bkg_peak_ch,
        &gain_table,
        hv_table.as_ref(),
        args.layer1_only,
        args.single_only,
    )?;

    let mut shifted_spectrum_fits = None;
    let mut lcbk_shift_fraction = None;
    if let Some(pha) = &bkg_pha {
        let cmap = ChannelMap::from_caldb_for_day(&args.caldb, unit, mission_day)
            .with_context(|| {
                format!(
                    "failed to load channel map for {} (day {:.1}) from {}",
                    unit.detnam(),
                    mission_day,
                    args.caldb.display()
                )
            })?;
        let shifted = shift_background_spectrum(
            unit,
            gain_shift.delta_ch_30kev,
            gain_shift.secular_scale,
            mission_day,
            &cmap,
            &pha.rates,
        )?;
        let filter = EventFilter {
            single_only: args.single_only,
            anodes: if args.layer1_only {
                AnodeSelection::Layer1
            } else {
                AnodeSelection::AllMain
            },
            pha_min: args.pha_min,
            pha_max: args.pha_max,
            ..Default::default()
        };
        validate_pha_range(filter.pha_min, filter.pha_max)?;

        let nch = shifted.len().max(1);
        let ncf = (1024 / nch).max(1);
        let ch_lo = (args.pha_min as usize) / ncf;
        let ch_hi = (args.pha_max as usize) / ncf;
        let (_cdel, fshbk) = background_lightcurve_shift_fraction(
            &pha.rates,
            &shifted,
            &filter.anodes,
            ch_lo,
            ch_hi,
        );
        lcbk_shift_fraction = Some(fshbk);

        if let Some(out_path) = &args.output {
            let mut rates = Vec::with_capacity(shifted.len());
            let mut errs = Vec::with_capacity(shifted.len());
            for row in &shifted {
                let s: f64 = (1..=10u8)
                    .filter(|&a| filter.anodes.matches(a))
                    .map(|a| row[(a - 1) as usize])
                    .sum();
                rates.push(s);
                errs.push((s / pha.exposure.max(1.0)).max(0.0).sqrt());
            }
            let mut meta = OgipMeta {
                instrume: "LAXPC".into(),
                detnam: Some(unit.detnam().into()),
                object: Some("BACKGROUND".into()),
                exposure: pha.exposure,
                ..Default::default()
            };
            meta.tstart = (mission_day + 2096.0) * 86_400.0;
            meta.tstop = meta.tstart + pha.exposure;
            ensure_parent_dir(out_path)?;
            write_spectrum_fits(
                out_path,
                unit,
                &rates,
                &errs,
                pha.exposure,
                43e-6,
                &filter,
                &gain_shift.recommended_rmf,
                "NONE",
                &meta,
            )?;
            shifted_spectrum_fits = Some(out_path.clone());
        }
    }

    Ok(LaxpcBackshiftReport {
        unit: unit.detnam().to_string(),
        mission_day,
        gain_shift,
        lcbk_shift_fraction,
        shifted_spectrum_fits,
    })
}

fn run_laxpc_orbit(args: &LaxpcOrbitArgs) -> Result<LaxpcOrbitReport> {
    let mkf = load_merged_mkf_files(&args.mkf)?;
    let tstart = mkf.times().first().copied().unwrap_or(0.0);
    let tstop = mkf.times().last().copied().unwrap_or(0.0);
    let mut meta = read_mkf_meta(&args.mkf, LaxpcUnit::Lxp10);
    meta.detnam = None;
    meta.tstart = tstart;
    meta.tstop = tstop;
    ensure_parent_dir(&args.output)?;
    write_orbit_from_mkf(&args.output, &mkf, &meta)?;

    Ok(LaxpcOrbitReport {
        mkf_files: args.mkf.len(),
        output: args.output.clone(),
        tstart,
        tstop,
    })
}

fn parse_units(raw: &[String]) -> Result<Vec<LaxpcUnit>> {
    let mut out = Vec::new();
    for item in raw {
        let t = item.trim();
        if t.eq_ignore_ascii_case("all") {
            return Ok(LaxpcUnit::ALL.to_vec());
        }
        let u = LaxpcUnit::from_str(t).map_err(|e| anyhow::anyhow!("{e}"))?;
        if !out.contains(&u) {
            out.push(u);
        }
    }
    if out.is_empty() {
        bail!("--units list cannot be empty");
    }
    Ok(out)
}

fn parse_anodes(raw: &str) -> Result<AnodeSelection> {
    let t = raw.trim().to_ascii_lowercase();
    if t == "all" || t == "a1-a7" {
        return Ok(AnodeSelection::AllMain);
    }
    if t == "layer1" || t == "l1" {
        return Ok(AnodeSelection::Layer1);
    }
    if t == "layer2" || t == "l2" {
        return Ok(AnodeSelection::Layer2);
    }
    let mut list = Vec::new();
    for part in t.split(',') {
        let p = part.trim().trim_start_matches('a');
        let id: u8 = p
            .parse()
            .with_context(|| format!("invalid anode ID '{part}' in --anodes"))?;
        if !(1..=10).contains(&id) {
            bail!("anode ID must be in 1..=10, got {id}");
        }
        if !list.contains(&id) {
            list.push(id);
        }
    }
    if list.len() == 1 {
        Ok(AnodeSelection::SingleAnode(list[0]))
    } else {
        Ok(AnodeSelection::Custom(list))
    }
}

fn resolve_laxpc_unit(explicit: Option<&str>, ea_paths: &[PathBuf]) -> Result<LaxpcUnit> {
    if let Some(s) = explicit {
        return LaxpcUnit::from_str(s).map_err(|e| anyhow::anyhow!("{e}"));
    }
    for p in ea_paths {
        let s = p.to_string_lossy().to_ascii_lowercase();
        if s.contains("lxp1") || s.contains("lx10") {
            return Ok(LaxpcUnit::Lxp10);
        }
        if s.contains("lxp2") || s.contains("lx20") {
            return Ok(LaxpcUnit::Lxp20);
        }
        if s.contains("lxp3") || s.contains("lx30") {
            return Ok(LaxpcUnit::Lxp30);
        }
    }
    bail!("could not infer LAXPC unit from filename; pass --unit <10|20|30>")
}

fn load_user_gti(path: &Path) -> Result<GtiSet> {
    if let Ok(gti) = GtiSet::from_fits(path, None) {
        return Ok(gti);
    }
    let txt = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read GTI file {}", path.display()))?;
    let mut intervals = Vec::new();
    for (idx, line) in txt.lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split_whitespace().collect();
        if cols.len() < 2 {
            bail!("GTI line {}: expected 'START STOP', got '{}'", idx + 1, t);
        }
        let a: f64 = cols[0].parse().with_context(|| format!("GTI line {} start", idx + 1))?;
        let b: f64 = cols[1].parse().with_context(|| format!("GTI line {} stop", idx + 1))?;
        intervals.push(GtiInterval::new(a, b));
    }
    if intervals.is_empty() {
        bail!("no GTI intervals found in {}", path.display());
    }
    Ok(GtiSet::new(intervals))
}

fn read_mkf_meta(mkf_paths: &[PathBuf], unit: LaxpcUnit) -> OgipMeta {
    let mut meta = OgipMeta {
        instrume: "LAXPC".to_string(),
        detnam: Some(unit.detnam().to_string()),
        ..Default::default()
    };
    if let Some(first) = mkf_paths.first() {
        if let Ok(fits) = FitsFile::open(first) {
            for idx in [0usize, 1usize] {
                if let Ok(h) = fits.hdu(idx) {
                    let hdr = h.header();
                    if meta.object.is_none() {
                        meta.object = hdr.get_str("OBJECT").map(|s| s.trim().to_string());
                    }
                    if meta.obs_id.is_none() {
                        meta.obs_id = hdr.get_str("OBS_ID").map(|s| s.trim().to_string());
                    }
                    if meta.ra_obj.is_none() {
                        meta.ra_obj = hdr
                            .get_f64("RA_NOM")
                            .or_else(|| hdr.get_f64("RA_OBJ"))
                            .or_else(|| hdr.get_f64("RA_PNT"));
                    }
                    if meta.dec_obj.is_none() {
                        meta.dec_obj = hdr
                            .get_f64("DEC_NOM")
                            .or_else(|| hdr.get_f64("DEC_OBJ"))
                            .or_else(|| hdr.get_f64("DEC_PNT"));
                    }
                }
            }
        }
    }
    meta
}

fn validate_pha_range(pha_min: u16, pha_max: u16) -> Result<()> {
    if pha_max > 1023 || pha_min > pha_max {
        bail!(
            "invalid PHA channel range [{pha_min}, {pha_max}]: require 0 <= pha_min <= pha_max <= 1023"
        );
    }
    Ok(())
}

fn merge_diag(dst: &mut FrameDiagnostics, src: FrameDiagnostics) {
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

fn parse_energy_range(er: &str) -> Result<(f64, f64)> {
    let (elo, ehi) = er
        .split_once(':')
        .context("--energy-range must be formatted as EMIN:EMAX (e.g. 3.0:80.0)")?;
    let lo: f64 = elo.trim().parse().context("invalid EMIN in --energy-range")?;
    let hi: f64 = ehi.trim().parse().context("invalid EMAX in --energy-range")?;
    if !lo.is_finite() || !hi.is_finite() || hi <= lo {
        bail!("--energy-range EMAX ({hi}) must be > EMIN ({lo})");
    }
    Ok((lo, hi))
}

fn validate_channel_maps(caldb_dir: Option<&Path>, units: &[LaxpcUnit]) -> Result<()> {
    let dir = caldb_dir.context(
        "--energy-range requires --caldb <dir> containing lxp{10,20,30}.channel files",
    )?;
    for &u in units {
        let cmap_path = dir.join(format!("{}0.channel", u.prefix()));
        if !cmap_path.is_file() {
            bail!(
                "--energy-range requires channel map {}, which was not found",
                cmap_path.display()
            );
        }
    }
    Ok(())
}

fn interp_lon_deg(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    if xs.is_empty() || ys.is_empty() {
        return 0.0;
    }
    let n = xs.len().min(ys.len());
    let norm = |d: f64| d.rem_euclid(360.0);
    if n == 1 || x <= xs[0] {
        return norm(ys[0]);
    }
    if x >= xs[n - 1] {
        return norm(ys[n - 1]);
    }
    let idx = xs[..n].partition_point(|&v| v <= x);
    let i0 = idx.saturating_sub(1);
    let i1 = idx.min(n - 1);
    let dx = xs[i1] - xs[i0];
    if dx.abs() <= 0.0 {
        return norm(ys[i0]);
    }
    let frac = (x - xs[i0]) / dx;
    let y0 = norm(ys[i0]);
    let mut y1 = norm(ys[i1]);
    if (y1 - y0).abs() > 180.0 {
        if y1 > y0 {
            y1 -= 360.0;
        } else {
            y1 += 360.0;
        }
    }
    norm(y0 + frac * (y1 - y0))
}

fn interp_linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    if xs.is_empty() || ys.is_empty() {
        return 0.0;
    }
    let n = xs.len().min(ys.len());
    if n == 1 || x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    let idx = xs[..n].partition_point(|&v| v <= x);
    let i0 = idx.saturating_sub(1);
    let i1 = idx.min(n - 1);
    let dx = xs[i1] - xs[i0];
    if dx.abs() <= 0.0 {
        ys[i0]
    } else {
        ys[i0] + (x - xs[i0]) / dx * (ys[i1] - ys[i0])
    }
}

fn ensure_parent_dir(path: &Path) -> Result<()> {
    if let Some(p) = path.parent() {
        if !p.as_os_str().is_empty() {
            std::fs::create_dir_all(p)?;
        }
    }
    Ok(())
}

fn print_laxpc_run_human(res: &MultiUnitPipelineResult) {
    println!(
        "LAXPC Level-1 -> Level-2 Pipeline Complete ({} unit(s)):",
        res.unit_results.len()
    );
    for u in &res.unit_results {
        println!(
            "  [{:<4}] orbits={}, valid_evts={}, gti_evts={}, net_exp={:.2} s, mean_rate={:.2} cps",
            u.unit.detnam(),
            u.num_orbits,
            u.total_valid_events,
            u.gti_filtered_events,
            u.gti.net_exposure_sec,
            u.lightcurve_mean_rate_cps
        );
        if let Some(ref rmf) = u.recommended_rmf {
            println!("         RMF recommendation: {rmf}");
        }
        println!(
            "         Outputs: spec={}, lc={}, gti={}{}",
            u.spectrum_fits_path.display(),
            u.lightcurve_fits_path.display(),
            u.gti_fits_path.display(),
            u.event_fits_path
                .as_ref()
                .map(|p| format!(", evn={}", p.display()))
                .unwrap_or_default()
        );
    }
    if let Some(ref orb) = res.orbit_fits_path {
        println!("  Orbit FITS       : {}", orb.display());
    }
    if let Some(ref fa) = res.format_a_event_fits_path {
        println!(
            "  Format A Events  : {} ({} events)",
            fa.display(),
            res.format_a_total_events
        );
    }
}
