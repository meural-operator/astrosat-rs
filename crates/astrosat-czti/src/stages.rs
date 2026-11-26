//! Standalone file-level CZTI stage entrypoints for `astrosat-cli` subcommands:
//! `bunchclean`, `pha2energy`, `gtigen`, `gaas`, `datasel`, `pixclean`, `noiseclean`,
//! `evtencor`, `evtclean`, `dpigen`, `image`, `bindata`.

use std::fs;
use std::path::{Path, PathBuf};

use astrosat_core::MkfTable;
use astrosat_fits::FitsFile;
use rayon::prelude::*;

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
    extract_mkf_rpy_rad, read_aux_tables, read_badpix_fits, read_czti_gtis, read_l2_event_file,
    read_mkf_bindata_attitude, read_response_file, read_spectrum_file, read_temp_extension,
    write_aspect_file, write_badpix_file, write_bunch_file, write_dph_file, write_dpi_file,
    write_gti_file, write_l2_event_file, write_lightcurve_file, write_livetime_file,
    write_response_file, write_sky_image_file, write_spectrum_file, CztiCaldb, CztiObsMeta,
    EventFileWriteOpts,
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
use crate::types::{
    detid_pixid_to_detxy, BadpixTable, BunchTable, GtiType, QuadAspect, QuadEvents, QuadLivetime,
    NUM_QUADS, PIXELS_PER_QUAD, QUAD_DIM,
};

/// Read a `*_bunch.fits` file (`Q0..Q3`).
pub fn read_bunch_file(path: impl AsRef<Path>) -> Result<[BunchTable; NUM_QUADS]> {
    let fits = FitsFile::open(path)?;
    let mut out = Vec::with_capacity(NUM_QUADS);
    for q in 0..NUM_QUADS {
        let mut b = BunchTable::new(q as u8);
        if let Ok(hdu) = fits.hdu_by_name(&format!("Q{q}")) {
            if let Ok(bt) = hdu.as_bintable() {
                b.time = bt.read_col_f64("Time").unwrap_or_default();
                b.time_dfs = bt.read_col_u8("Time_dfs").unwrap_or_default();
                b.time_dsl = bt.read_col_u8("Time_dsl").unwrap_or_default();
                b.num_event = bt.read_col_u8("NumEvent").unwrap_or_default();
                b.detid1 = bt.read_col_u8("DetId1").unwrap_or_default();
                b.detid2 = bt.read_col_u8("DetId2").unwrap_or_default();
                b.detid3 = bt.read_col_u8("DetId3").unwrap_or_default();
                b.detid4 = bt.read_col_u8("DetId4").unwrap_or_default();
                b.row_num = bt.read_col_i32("RowNum").unwrap_or_default();
                b.detid_fevt = bt.read_col_u8("DetId_fevt").unwrap_or_default();
                b.pixid_fevt = bt.read_col_u8("PixId_fevt").unwrap_or_default();
                b.detid_sevt = bt.read_col_u8("DetId_sevt").unwrap_or_default();
                b.pixid_sevt = bt.read_col_u8("PixId_sevt").unwrap_or_default();
                b.detid_levt = bt.read_col_u8("DetId_levt").unwrap_or_default();
                b.pixid_levt = bt.read_col_u8("PixId_levt").unwrap_or_default();
            }
        }
        out.push(b);
    }
    let mut it = out.into_iter();
    Ok([
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
    ])
}

/// Read a `*_livetime.fits` file (`Q0..Q3`).
pub fn read_livetime_file(path: impl AsRef<Path>) -> Result<[QuadLivetime; NUM_QUADS]> {
    let fits = FitsFile::open(path)?;
    let mut out = Vec::with_capacity(NUM_QUADS);
    for q in 0..NUM_QUADS {
        let hdu = fits.hdu_by_name(&format!("Q{q}")).or_else(|_| fits.hdu(q + 1))?;
        let bt = hdu.as_bintable()?;
        let time = bt.read_col_f64("TIME")?;
        let fracexp = bt.read_col_f64("FRACEXP")?;
        let binsize = if time.len() >= 2 {
            (time[1] - time[0]).max(1e-6)
        } else {
            hdu.header().get_f64("LV_BINSIZE").unwrap_or(1.0)
        };
        out.push(QuadLivetime {
            quad_id: q as u8,
            binsize,
            time,
            fracexp,
        });
    }
    let mut it = out.into_iter();
    Ok([
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
    ])
}

/// Read an aspect file (`*.aspect_Q{q}`).
pub fn read_aspect_file(path: impl AsRef<Path>, quad_id: u8) -> Result<QuadAspect> {
    let fits = FitsFile::open(path)?;
    let hdu = fits.hdu_by_name("ASPECT").or_else(|_| fits.hdu(1))?;
    let hdr = hdu.header();
    let bt = hdu.as_bintable()?;
    Ok(QuadAspect {
        quad_id,
        teldef_name: hdr.get_str("TELDEF").unwrap_or_default().trim().to_string(),
        time: bt.read_col_f64("TIME").unwrap_or_default(),
        nx: bt.read_col_f32("Nx").unwrap_or_default(),
        ny: bt.read_col_f32("Ny").unwrap_or_default(),
        nz: bt.read_col_f32("Nz").unwrap_or_default(),
        nxt: bt.read_col_f32("Nxt").unwrap_or_default(),
        nyt: bt.read_col_f32("Nyt").unwrap_or_default(),
        nzt: bt.read_col_f32("Nzt").unwrap_or_default(),
        ra_deg: hdr.get_f64("RA").unwrap_or(0.0),
        dec_deg: hdr.get_f64("DEC").unwrap_or(0.0),
        twist_deg: hdr.get_f64("TWIST").unwrap_or(0.0),
    })
}

/// Standalone Stage 1b: `cztbunchclean`.
pub fn run_bunchclean(
    in_evt: impl AsRef<Path>,
    in_bunch: Option<&Path>,
    out_evt: impl AsRef<Path>,
    out_livetime: impl AsRef<Path>,
    out_bunch: impl AsRef<Path>,
    cfg: &BunchCleanConfig,
) -> Result<[usize; NUM_QUADS]> {
    let (meta, quads, gtis, _, _) = read_l2_event_file(&in_evt)?;
    let bunches = match in_bunch {
        Some(p) if p.exists() => read_bunch_file(p)?,
        _ => [
            BunchTable::new(0),
            BunchTable::new(1),
            BunchTable::new(2),
            BunchTable::new(3),
        ],
    };

    let res: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| bunchclean_quad(&quads[q], &bunches[q], &gtis.quad[q], cfg))
        .collect();

    let counts = [
        res[0].events.len(),
        res[1].events.len(),
        res[2].events.len(),
        res[3].events.len(),
    ];
    let ncompat = res[0].ncompat;
    let out_q = [
        res[0].events.clone(),
        res[1].events.clone(),
        res[2].events.clone(),
        res[3].events.clone(),
    ];
    let out_b = [
        res[0].out_bunches.clone(),
        res[1].out_bunches.clone(),
        res[2].out_bunches.clone(),
        res[3].out_bunches.clone(),
    ];
    let out_lv = [
        res[0].livetime.clone(),
        res[1].livetime.clone(),
        res[2].livetime.clone(),
        res[3].livetime.clone(),
    ];
    let pexp = [
        res[0].pixel_exposure.clone(),
        res[1].pixel_exposure.clone(),
        res[2].pixel_exposure.clone(),
        res[3].pixel_exposure.clone(),
    ];

    let (aux_veto, aux_ssm, aux_temp) = read_aux_tables(&in_evt);
    write_l2_event_file(
        out_evt,
        &meta,
        &out_q,
        &aux_veto,
        &aux_ssm,
        &aux_temp,
        &gtis,
        &EventFileWriteOpts {
            ncompat: Some(ncompat),
            pixel_exposure: Some(pexp),
            ..Default::default()
        },
    )?;
    write_livetime_file(out_livetime, &meta, &out_lv)?;
    write_bunch_file(out_bunch, &meta, &out_b)?;
    Ok(counts)
}

/// Standalone Stage 1c: `cztpha2energy`.
pub fn run_pha2energy(
    in_evt: impl AsRef<Path>,
    temp_ss_file: Option<&Path>,
    caldb_dir: impl AsRef<Path>,
    out_evt: impl AsRef<Path>,
) -> Result<[usize; NUM_QUADS]> {
    let (meta, mut quads, gtis, pexp, _) = read_l2_event_file(&in_evt)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let temp = temp_ss_file
        .and_then(|p| read_temp_extension(p).ok())
        .or_else(|| read_temp_extension(&in_evt).ok())
        .unwrap_or_default();

    let ebounds = caldb.read_ebounds(meta.tstart)?;
    let gains: Vec<_> = (0..NUM_QUADS)
        .map(|q| caldb.read_gain_offset(q, meta.tstart))
        .collect::<Result<Vec<_>>>()?;
    let tseries: Vec<_> = (0..NUM_QUADS)
        .map(|q| QuadTempSeries::from_temp_table(&temp, q as u8))
        .collect();

    quads
        .par_iter_mut()
        .enumerate()
        .for_each(|(q, ev)| pha2energy_quad(ev, &tseries[q], &gains[q], &ebounds));

    let (aux_veto, aux_ssm, _) = read_aux_tables(&in_evt);
    let counts = [quads[0].len(), quads[1].len(), quads[2].len(), quads[3].len()];
    write_l2_event_file(
        out_evt,
        &meta,
        &quads,
        &aux_veto,
        &aux_ssm,
        &temp,
        &gtis,
        &EventFileWriteOpts {
            pixel_exposure: pexp,
            ..Default::default()
        },
    )?;
    Ok(counts)
}

/// Standalone Stage 2a: `cztgtigen`.
pub fn run_gtigen(
    in_evt: impl AsRef<Path>,
    mkf_file: impl AsRef<Path>,
    threshold_file: Option<&Path>,
    user_gti_file: Option<&Path>,
    out_gti: impl AsRef<Path>,
    cfg: &GtiGenConfig,
) -> Result<[f64; NUM_QUADS]> {
    let meta = CztiObsMeta::from_fits(&in_evt)?;
    let (m0_gtis, m0_bti) = read_czti_gtis(&in_evt)?;
    let mkf = MkfTable::from_fits(mkf_file)?;
    let rules = if let Some(tf) = threshold_file {
        parse_mkf_threshold_text(&fs::read_to_string(tf)?)
    } else {
        cfg.thresholds.clone()
    };

    let mut good = vec![true; mkf.len()];
    for rule in &rules {
        if !rule.enabled {
            continue;
        }
        if let Some(col) = mkf.column(&rule.column) {
            for (i, &v) in col.iter().enumerate() {
                if !v.is_finite()
                    || rule.min.is_some_and(|lo| v < lo)
                    || rule.max.is_some_and(|hi| v > hi)
                {
                    good[i] = false;
                }
            }
        }
    }
    let mkf_gti = intervals_from_mkf_mask(mkf.times(), &good);
    let user_gti = match user_gti_file {
        Some(p) => Some(parse_user_gti_text(&fs::read_to_string(p)?)?),
        None => None,
    };

    let out = generate_czti_gtis(
        &mkf_gti,
        &m0_gtis.common,
        &m0_gtis.quad,
        if m0_bti.is_empty() { None } else { Some(&m0_bti) },
        user_gti.as_deref(),
        cfg.min_exposure,
    )?;
    write_gti_file(out_gti, &meta, &out)?;

    let mut exps = [0.0f64; NUM_QUADS];
    for q in 0..NUM_QUADS {
        exps[q] = out.quad[q].iter().map(|&(a, b)| (b - a).max(0.0)).sum();
    }
    Ok(exps)
}

/// Standalone Stage 2b: `cztgaas`.
pub fn run_gaas(
    in_evt: impl AsRef<Path>,
    mkf_file: impl AsRef<Path>,
    caldb_dir: impl AsRef<Path>,
    out_aspect_prefix: impl AsRef<Path>,
) -> Result<[QuadAspect; NUM_QUADS]> {
    let (meta, quads, gtis, _, _) = read_l2_event_file(&in_evt)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let mkf = MkfTable::from_fits(mkf_file)?;
    let rpy = extract_mkf_rpy_rad(&mkf)?;

    let mut aspects = Vec::with_capacity(NUM_QUADS);
    let prefix_str = out_aspect_prefix.as_ref().to_string_lossy().to_string();
    for q in 0..NUM_QUADS {
        let (tname, tmat) = caldb.read_teldef(q, meta.tstart)?;
        let t0 = if quads[q].tstart > 0.0 { quads[q].tstart } else { meta.tstart };
        let t1 = if quads[q].tstop > 0.0 { quads[q].tstop } else { meta.tstop };
        let asp = compute_quad_aspect(
            q as u8,
            tname,
            &tmat,
            &rpy,
            t0,
            t1,
            if gtis.quad[q].is_empty() { None } else { Some(&gtis.quad[q]) },
        );
        write_aspect_file(format!("{prefix_str}_Q{q}"), &meta, &asp)?;
        aspects.push(asp);
    }
    let mut it = aspects.into_iter();
    Ok([
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
    ])
}

/// Standalone Stage 2c: `cztdatasel`.
pub fn run_datasel(
    in_evt: impl AsRef<Path>,
    gti_file: impl AsRef<Path>,
    gti_type: GtiType,
    out_evt: impl AsRef<Path>,
) -> Result<[usize; NUM_QUADS]> {
    let (meta, quads, _, pexp, gcor) = read_l2_event_file(&in_evt)?;
    let (gtis, _) = read_czti_gtis(&gti_file)?;
    let filtered: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| {
            let ivs = match gti_type {
                GtiType::Common => &gtis.common,
                GtiType::Quad => &gtis.quad[q],
            };
            datasel_quad(&quads[q], ivs)
        })
        .collect();
    let counts = [
        filtered[0].len(),
        filtered[1].len(),
        filtered[2].len(),
        filtered[3].len(),
    ];
    let arr = [
        filtered[0].clone(),
        filtered[1].clone(),
        filtered[2].clone(),
        filtered[3].clone(),
    ];
    let (aux_veto, aux_ssm, aux_temp) = read_aux_tables(&in_evt);
    write_l2_event_file(
        out_evt,
        &meta,
        &arr,
        &aux_veto,
        &aux_ssm,
        &aux_temp,
        &gtis,
        &EventFileWriteOpts {
            gti_type: Some(gti_type),
            gaincor_applied: gcor,
            pixel_exposure: pexp,
            ..Default::default()
        },
    )?;
    Ok(counts)
}

/// Standalone Stage 2d (Default): `cztpixclean`.
#[allow(clippy::too_many_arguments)]
pub fn run_pixclean(
    in_evt: impl AsRef<Path>,
    in_livetime: impl AsRef<Path>,
    caldb_dir: impl AsRef<Path>,
    out_evt: impl AsRef<Path>,
    out_livetime: impl AsRef<Path>,
    out_badpix: impl AsRef<Path>,
    gti_type: GtiType,
    cfg: &PixCleanConfig,
) -> Result<[usize; NUM_QUADS]> {
    let (meta, quads, gtis, pexp_opt, gcor) = read_l2_event_file(&in_evt)?;
    let in_lv = read_livetime_file(in_livetime)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let caldb_bp = caldb.read_badpix(meta.tstart)?;
    let in_pexp = pexp_opt.unwrap_or_else(|| {
        [
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
        ]
    });

    let res: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| {
            let ivs = match gti_type {
                GtiType::Common => &gtis.common,
                GtiType::Quad => &gtis.quad[q],
            };
            pixclean_quad(&quads[q], &in_lv[q], ivs, &caldb_bp.flags[q], &in_pexp[q], cfg)
        })
        .collect();

    let counts = [
        res[0].single_events.len(),
        res[1].single_events.len(),
        res[2].single_events.len(),
        res[3].single_events.len(),
    ];
    let out_q = [
        res[0].single_events.clone(),
        res[1].single_events.clone(),
        res[2].single_events.clone(),
        res[3].single_events.clone(),
    ];
    let out_lv = [
        res[0].livetime.clone(),
        res[1].livetime.clone(),
        res[2].livetime.clone(),
        res[3].livetime.clone(),
    ];
    let mut bp = BadpixTable::default();
    let mut out_pexp: [Vec<f64>; NUM_QUADS] = [
        vec![0.0f64; PIXELS_PER_QUAD],
        vec![0.0f64; PIXELS_PER_QUAD],
        vec![0.0f64; PIXELS_PER_QUAD],
        vec![0.0f64; PIXELS_PER_QUAD],
    ];
    for q in 0..NUM_QUADS {
        bp.flags[q] = res[q].badpix_flags.clone();
        out_pexp[q] = res[q].pixel_exposure.clone();
    }

    let (aux_veto, aux_ssm, aux_temp) = read_aux_tables(&in_evt);
    write_l2_event_file(
        out_evt,
        &meta,
        &out_q,
        &aux_veto,
        &aux_ssm,
        &aux_temp,
        &gtis,
        &EventFileWriteOpts {
            gti_type: Some(gti_type),
            gaincor_applied: gcor,
            pixel_exposure: Some(out_pexp),
            ..Default::default()
        },
    )?;
    write_livetime_file(out_livetime, &meta, &out_lv)?;
    write_badpix_file(out_badpix, &meta, &bp)?;
    Ok(counts)
}

/// Standalone Stage 2d (Generalized v3.0): combined `noiseclean` (`noisypix -> superbunch ->
/// heavybunch -> flickpix -> eventsep`).
#[allow(clippy::too_many_arguments)]
pub fn run_noiseclean(
    in_evt: impl AsRef<Path>,
    in_bunch: impl AsRef<Path>,
    in_livetime: impl AsRef<Path>,
    caldb_dir: impl AsRef<Path>,
    out_evt: impl AsRef<Path>,
    out_livetime: impl AsRef<Path>,
    out_badpix: impl AsRef<Path>,
    cfg: &NoiseCleanConfig,
) -> Result<[usize; NUM_QUADS]> {
    let (meta, quads, mut gtis, pexp_opt, gcor) = read_l2_event_file(&in_evt)?;
    let (aux_veto, aux_ssm, aux_temp) = read_aux_tables(&in_evt);
    let bunches = read_bunch_file(in_bunch)?;
    let in_lv = read_livetime_file(in_livetime)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let caldb_bp = caldb.read_badpix(meta.tstart)?;
    let gaincor: Vec<_> = (0..NUM_QUADS)
        .map(|q| caldb.read_gaincor_lld_uld(q, meta.tstart))
        .collect::<Result<Vec<_>>>()?;
    let in_pexp = pexp_opt.unwrap_or_else(|| {
        [
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
        ]
    });

    let res: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| {
            let np = noisypixclean_quad(&quads[q], &caldb_bp.flags[q], &in_pexp[q], cfg);
            let sb = superbunchclean_quad(&np.events, &bunches[q], &in_lv[q], cfg);
            let (hb_ev, hb_pexp) = heavybunchclean_quad_with_exposure(
                &sb.events,
                &bunches[q],
                &gaincor[q].lld,
                &np.pixel_exposure,
                cfg,
            );
            let fp = flickpixclean_quad_with_gti(
                &hb_ev,
                &np.badpix_flags,
                &hb_pexp,
                Some(&gtis.quad[q]),
                cfg,
            );
            let (mut singles, doubles) = eventsep_quad(&fp.events, cfg.double_time_sec);
            if !sb.bti_intervals.is_empty() {
                let updated_gti = subtract_intervals(&gtis.quad[q], &sb.bti_intervals);
                singles.exposure = updated_gti.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
            }
            (
                singles,
                doubles,
                sb.livetime,
                fp.badpix_flags,
                fp.pixel_exposure,
                sb.bti_intervals,
            )
        })
        .collect();

    let counts = [res[0].0.len(), res[1].0.len(), res[2].0.len(), res[3].0.len()];
    let out_q = [
        res[0].0.clone(),
        res[1].0.clone(),
        res[2].0.clone(),
        res[3].0.clone(),
    ];
    let out_dbl = [
        res[0].1.clone(),
        res[1].1.clone(),
        res[2].1.clone(),
        res[3].1.clone(),
    ];
    let out_lv = [
        res[0].2.clone(),
        res[1].2.clone(),
        res[2].2.clone(),
        res[3].2.clone(),
    ];
    let mut bp = BadpixTable::default();
    let mut out_pexp: [Vec<f64>; NUM_QUADS] = [
        vec![0.0f64; PIXELS_PER_QUAD],
        vec![0.0f64; PIXELS_PER_QUAD],
        vec![0.0f64; PIXELS_PER_QUAD],
        vec![0.0f64; PIXELS_PER_QUAD],
    ];
    for q in 0..NUM_QUADS {
        bp.flags[q] = res[q].3.clone();
        out_pexp[q] = res[q].4.clone();
        if !res[q].5.is_empty() {
            gtis.quad[q] = subtract_intervals(&gtis.quad[q], &res[q].5);
        }
    }

    write_l2_event_file(
        &out_evt,
        &meta,
        &out_q,
        &aux_veto,
        &aux_ssm,
        &aux_temp,
        &gtis,
        &EventFileWriteOpts {
            gti_type: Some(GtiType::Quad),
            ncompat: Some(1),
            gaincor_applied: gcor,
            pixel_exposure: Some(out_pexp.clone()),
            ..Default::default()
        },
    )?;
    if cfg.write_double_evt {
        let dbl_path = out_evt.as_ref().with_extension("double.evt");
        write_l2_event_file(
            &dbl_path,
            &meta,
            &out_dbl,
            &aux_veto,
            &aux_ssm,
            &aux_temp,
            &gtis,
            &EventFileWriteOpts {
                gti_type: Some(GtiType::Quad),
                ncompat: Some(1),
                gaincor_applied: gcor,
                pixel_exposure: Some(out_pexp),
                ..Default::default()
            },
        )?;
    }
    write_livetime_file(out_livetime, &meta, &out_lv)?;
    write_badpix_file(out_badpix, &meta, &bp)?;
    Ok(counts)
}

/// Standalone Stage 2e: `cztevtencor` (idempotent `GAINCOR='YES'` guard).
pub fn run_evtencor(
    in_evt: impl AsRef<Path>,
    caldb_dir: impl AsRef<Path>,
    out_evt: impl AsRef<Path>,
    cfg: &EvtEnCorConfig,
) -> Result<[usize; NUM_QUADS]> {
    let (meta, quads, gtis, pexp, already_corrected) = read_l2_event_file(&in_evt)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let gcor_caldb: Vec<_> = (0..NUM_QUADS)
        .map(|q| caldb.read_gaincor_lld_uld(q, meta.tstart))
        .collect::<Result<Vec<_>>>()?;

    let res: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| evtencor_quad(&quads[q], &gcor_caldb[q], meta.tstart, already_corrected, cfg))
        .collect();

    let counts = [
        res[0].events.len(),
        res[1].events.len(),
        res[2].events.len(),
        res[3].events.len(),
    ];
    let out_q = [
        res[0].events.clone(),
        res[1].events.clone(),
        res[2].events.clone(),
        res[3].events.clone(),
    ];
    let (aux_veto, aux_ssm, aux_temp) = read_aux_tables(&in_evt);
    write_l2_event_file(
        out_evt,
        &meta,
        &out_q,
        &aux_veto,
        &aux_ssm,
        &aux_temp,
        &gtis,
        &EventFileWriteOpts {
            gaincor_applied: true,
            lld_cut: cfg.apply_lld,
            uld_cut: cfg.apply_uld,
            pixel_exposure: pexp,
            ..Default::default()
        },
    )?;
    Ok(counts)
}

/// Standalone Stage 2f: `cztevtclean`.
pub fn run_evtclean(
    in_evt: impl AsRef<Path>,
    out_evt: impl AsRef<Path>,
    cfg: &EvtCleanConfig,
) -> Result<[usize; NUM_QUADS]> {
    let (meta, quads, gtis, pexp, gcor) = read_l2_event_file(&in_evt)?;
    let cleaned: Vec<_> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| evtclean_quad(&quads[q], cfg))
        .collect();
    let counts = [
        cleaned[0].len(),
        cleaned[1].len(),
        cleaned[2].len(),
        cleaned[3].len(),
    ];
    let out_q = [
        cleaned[0].clone(),
        cleaned[1].clone(),
        cleaned[2].clone(),
        cleaned[3].clone(),
    ];
    let (aux_veto, aux_ssm, aux_temp) = read_aux_tables(&in_evt);
    write_l2_event_file(
        out_evt,
        &meta,
        &out_q,
        &aux_veto,
        &aux_ssm,
        &aux_temp,
        &gtis,
        &EventFileWriteOpts {
            gaincor_applied: gcor,
            pixel_exposure: pexp,
            ..Default::default()
        },
    )?;
    Ok(counts)
}

/// Standalone Stage 3a: `cztdpigen`.
pub fn run_dpigen(
    in_evt: impl AsRef<Path>,
    badpix_file: impl AsRef<Path>,
    caldb_dir: impl AsRef<Path>,
    out_dph: impl AsRef<Path>,
    out_dpi: impl AsRef<Path>,
    cfg: &DpiGenConfig,
) -> Result<()> {
    let (meta, quads, _, _, _) = read_l2_event_file(&in_evt)?;
    let badpix = read_badpix_fits(badpix_file)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let eff_areas = [
        caldb.read_eff_area(0, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(1, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(2, meta.tstart).unwrap_or_default(),
        caldb.read_eff_area(3, meta.tstart).unwrap_or_default(),
    ];
    let prod = generate_dph_dpi(&quads, &badpix, &eff_areas, cfg);
    write_dph_file(out_dph, &meta, GtiType::Quad, &prod.quad_dph, &prod.full_dph)?;
    write_dpi_file(out_dpi, &meta, GtiType::Quad, &prod.quad_dpi, &prod.full_dpi)?;
    Ok(())
}

/// Standalone Stage 3b: `cztimage`.
pub fn run_image(
    in_dpi: impl AsRef<Path>,
    aspect_files: &[PathBuf; NUM_QUADS],
    caldb_dir: impl AsRef<Path>,
    out_img_prefix: impl AsRef<Path>,
    cfg: &ImageConfig,
) -> Result<Vec<PathBuf>> {
    let meta = CztiObsMeta::from_fits(&in_dpi)?;
    let fits = FitsFile::open(&in_dpi)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let prefix_str = out_img_prefix.as_ref().to_string_lossy().to_string();
    let mut written = Vec::new();

    let shadowlib_all = if cfg.algorithm == ImageAlgorithm::BalCC {
        caldb.read_shadowlib_all(meta.tstart).ok()
    } else {
        None
    };

    let mut sky_imgs = Vec::new();
    let mut first_aspect = None;
    for &q in &cfg.quads {
        if q >= NUM_QUADS {
            continue;
        }
        let hdu = fits.hdu_by_name(&format!("Q{q}")).or_else(|_| fits.hdu(q + 1))?;
        let flat = hdu.as_image()?.read_f32()?;
        let mut dpi = vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM];
        for y in 0..QUAD_DIM {
            for x in 0..QUAD_DIM {
                if y * QUAD_DIM + x < flat.len() {
                    dpi[y][x] = flat[y * QUAD_DIM + x];
                }
            }
        }
        let aspect = read_aspect_file(&aspect_files[q], q as u8)?;
        if first_aspect.is_none() {
            first_aspect = Some(aspect.clone());
        }
        let sky = match cfg.algorithm {
            ImageAlgorithm::BalCC => {
                if let Some((grid_dim, ref all_shadows)) = shadowlib_all {
                    let refs: Vec<&[f32]> = all_shadows[q].iter().map(|v| v.as_slice()).collect();
                    balanced_cross_correlation(&dpi, &refs, grid_dim, &aspect)
                } else {
                    let mask = caldb
                        .read_mask_pattern(q, meta.tstart)
                        .unwrap_or_else(|_| vec![vec![0.5; QUAD_DIM]; QUAD_DIM]);
                    balanced_cross_correlation_from_mask(&dpi, &mask, 15, &aspect)
                }
            }
            ImageAlgorithm::Fft => {
                let mask = caldb
                    .read_mask_pattern(q, meta.tstart)
                    .unwrap_or_else(|_| vec![vec![0.5; QUAD_DIM]; QUAD_DIM]);
                fft_cross_correlation(&dpi, &mask, cfg.oversampling_factor, &aspect)
            }
        };
        let out_path = PathBuf::from(format!("{prefix_str}_Q{q}"));
        write_sky_image_file(&out_path, &format!("Q{q}"), &meta, &sky)?;
        written.push(out_path);
        sky_imgs.push(sky);
    }
    if sky_imgs.len() > 1 {
        if let Some(ref asp0) = first_aspect {
            let combined = combine_sky_images(&sky_imgs, asp0);
            let comb_path = PathBuf::from(format!("{prefix_str}_all"));
            write_sky_image_file(&comb_path, "IMAGE", &meta, &combined)?;
            written.push(comb_path);
        }
    }
    Ok(written)
}

/// Standalone Stage 3c: `cztbindata`.
#[allow(clippy::too_many_arguments)]
pub fn run_bindata(
    in_evt: impl AsRef<Path>,
    mkf_file: impl AsRef<Path>,
    badpix_file: impl AsRef<Path>,
    livetime_file: impl AsRef<Path>,
    caldb_dir: impl AsRef<Path>,
    out_prefix: impl AsRef<Path>,
    cfg: &BinDataConfig,
) -> Result<Vec<PathBuf>> {
    let (meta, quads, gtis, pexp_opt, already_gcor) = read_l2_event_file(&in_evt)?;
    let (aux_veto, aux_ssm, aux_temp) = read_aux_tables(&in_evt);
    let badpix = read_badpix_fits(badpix_file)?;
    let livetimes = read_livetime_file(livetime_file)?;
    let caldb = CztiCaldb::open(caldb_dir)?;
    let ebounds = caldb.read_ebounds(meta.tstart)?;
    let pexp = pexp_opt.unwrap_or_else(|| {
        [
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
            vec![1.0; PIXELS_PER_QUAD],
        ]
    });

    let (ra_pnt, dec_pnt, att_t, roll_ra, roll_dec, roll_rot) = read_mkf_bindata_attitude(&mkf_file)?;
    let (ra_use, dec_use) = cfg.resolve_source_radec(ra_pnt, dec_pnt);
    let att_tx_ty: Vec<(f64, f64)> = (0..att_t.len())
        .into_par_iter()
        .map(|i| radec_to_czti_tx_ty(roll_ra[i], roll_dec[i], roll_rot[i], ra_use, dec_use))
        .collect();

    let prefix_str = out_prefix.as_ref().to_string_lossy().to_string();
    let mut written = Vec::new();
    let mut weighted_quads: [QuadEvents; NUM_QUADS] = quads.clone();
    let mut has_weighted = false;
    let mut common_specs: Vec<CztiSpectrum> = Vec::new();
    let mut common_lcs: Vec<CztiLightCurve> = Vec::new();
    let mut common_rsps: Vec<CztiResponse> = Vec::new();

    let shared_pixresp = if cfg.gen_resp
        && matches!(cfg.output_type, BinDataOutputType::Spec | BinDataOutputType::Both)
    {
        caldb.read_pixresp(meta.tstart).ok()
    } else {
        None
    };

    for &q in &cfg.quads {
        if q >= NUM_QUADS {
            continue;
        }
        let qgti = match cfg.gti_type {
            GtiType::Common => &gtis.common,
            GtiType::Quad => &gtis.quad[q],
        };
        let mut sum_tx = 0.0;
        let mut sum_ty = 0.0;
        let mut gtisel_tx_ty: Vec<(f64, f64)> = Vec::new();
        for (i, &t) in att_t.iter().enumerate() {
            if qgti.is_empty() || qgti.iter().any(|&(a, b)| t >= a && t <= b) {
                let (tx, ty) = att_tx_ty[i];
                if cfg.mask_weight || cfg.gen_resp {
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
        } else {
            (0.0, 0.0)
        };

        let mut gaincor = caldb.read_gaincor_lld_uld(q, meta.tstart)?;
        apply_additional_lld(q, &mut gaincor.lld);
        let mut lld_detmap = [[0i16; QUAD_DIM]; QUAD_DIM];
        let mut uld_detmap = [[511i16; QUAD_DIM]; QUAD_DIM];
        let mut pix_exp_detxy = vec![0.0f64; PIXELS_PER_QUAD];
        for det in 0..16u8 {
            for pix in 0..=255u8 {
                let row = (det as usize) * 256 + (pix as usize);
                let (dx, dy) = detid_pixid_to_detxy(q as u8, det, pix);
                lld_detmap[dx as usize][dy as usize] = gaincor.lld[row];
                uld_detmap[dx as usize][dy as usize] = gaincor.uld[row];
                pix_exp_detxy[(dx as usize) * QUAD_DIM + (dy as usize)] = pexp[q][row];
            }
        }

        let mut tx_sh = DETWISE_TX_SHIFT[q];
        let mut ty_sh = DETWISE_TY_SHIFT[q];
        for m in 0..16 {
            tx_sh[m] += cfg.delta_tx;
            ty_sh[m] += cfg.delta_ty;
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

        let renorm = if cfg.mask_weight {
            let bkgscal = caldb
                .read_bkgscal(q, meta.tstart)
                .unwrap_or_else(|_| vec![1.0; PIXELS_PER_QUAD * 512]);
            Some(calculate_renormalized_weights(
                q,
                &open_frac,
                &raw_w,
                &badpix,
                cfg.badpix_threshold,
                &bkgscal,
                &lld_detmap,
                &uld_detmap,
                &pix_exp_detxy,
                cfg.min_active_pix_per_module,
            ))
        } else {
            None
        };

        let syserr = caldb.read_syserr(q, meta.tstart).unwrap_or([0.0f32; 512]);
        let rsp_path = if cfg.gti_type == GtiType::Common {
            PathBuf::from(format!("{prefix_str}.rsp"))
        } else {
            PathBuf::from(format!("{prefix_str}_Q{q}.rsp"))
        };
        let rsp_name = if cfg.gen_resp {
            rsp_path
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "NONE".to_string())
        } else {
            "NONE".to_string()
        };

        let (spec, lc, wevt) = bindata_quad(
            &quads[q],
            &badpix,
            &livetimes[q],
            renorm.as_ref(),
            &gaincor,
            already_gcor,
            &syserr,
            meta.tstart,
            meta.tstop,
            avg_tx,
            avg_ty,
            rsp_name,
            cfg,
        );

        let quad_rsp = if cfg.gen_resp
            && matches!(cfg.output_type, BinDataOutputType::Spec | BinDataOutputType::Both)
        {
            if let Some(ref pixresp) = shared_pixresp {
                if let Ok(mut rsp_caldb) =
                    caldb.read_rsp_caldb_with_pixresp(q, meta.tstart, pixresp)
                {
                    apply_additional_lld(q, &mut rsp_caldb.lld);
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
                        &badpix,
                        cfg.badpix_threshold,
                        &avg_of_5x5,
                        renorm.as_ref(),
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

        if cfg.gti_type == GtiType::Common {
            common_specs.push(spec);
            common_lcs.push(lc);
            if let Some(r) = quad_rsp {
                common_rsps.push(r);
            }
        } else {
            if matches!(cfg.output_type, BinDataOutputType::Spec | BinDataOutputType::Both) {
                let pha_path = PathBuf::from(format!("{prefix_str}_Q{q}.pha"));
                write_spectrum_file(&pha_path, &meta, &spec)?;
                written.push(pha_path);

                if let Some(rsp) = quad_rsp {
                    write_response_file(&rsp_path, &meta, &rsp)?;
                    written.push(rsp_path);
                }
            }
            if matches!(cfg.output_type, BinDataOutputType::Lc | BinDataOutputType::Both) {
                let lc_path = PathBuf::from(format!("{prefix_str}_Q{q}.lc"));
                write_lightcurve_file(&lc_path, &meta, &lc)?;
                written.push(lc_path);
            }
        }
        if let Some(we) = wevt {
            weighted_quads[q] = we;
            has_weighted = true;
        }
    }

    if cfg.gti_type == GtiType::Common {
        let rsp_path = PathBuf::from(format!("{prefix_str}.rsp"));
        let rsp_name = if cfg.gen_resp && !common_rsps.is_empty() {
            rsp_path
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "NONE".to_string())
        } else {
            "NONE".to_string()
        };
        if matches!(cfg.output_type, BinDataOutputType::Spec | BinDataOutputType::Both)
            && !common_specs.is_empty()
        {
            let comb_spec = combine_spectra(&common_specs, false, rsp_name);
            let pha_path = PathBuf::from(format!("{prefix_str}.pha"));
            write_spectrum_file(&pha_path, &meta, &comb_spec)?;
            written.push(pha_path);

            if !common_rsps.is_empty() {
                let exps: Vec<f64> = common_specs.iter().map(|s| s.exposure).collect();
                let comb_rsp = combine_responses(&common_rsps, &exps, false);
                write_response_file(&rsp_path, &meta, &comb_rsp)?;
                written.push(rsp_path);
            }
        }
        if matches!(cfg.output_type, BinDataOutputType::Lc | BinDataOutputType::Both)
            && !common_lcs.is_empty()
        {
            let comb_lc = combine_lightcurves(&common_lcs);
            let lc_path = PathBuf::from(format!("{prefix_str}.lc"));
            write_lightcurve_file(&lc_path, &meta, &comb_lc)?;
            written.push(lc_path);
        }
    }

    if has_weighted {
        let wp = PathBuf::from(format!("{prefix_str}_weight.evt"));
        write_l2_event_file(
            &wp,
            &meta,
            &weighted_quads,
            &aux_veto,
            &aux_ssm,
            &aux_temp,
            &gtis,
            &EventFileWriteOpts {
                gti_type: Some(cfg.gti_type),
                gaincor_applied: true,
                pixel_exposure: Some(pexp),
                ..Default::default()
            },
        )?;
        written.push(wp);
    }

    Ok(written)
}

/// Standalone `cztaddspec` — co-add multiple CZTI OGIP `.pha` spectra and (optionally) exposure-weighted
/// `.rsp` response matrices (`cztaddspec.cpp`).
pub fn run_addspec(
    pha_files: &[PathBuf],
    rsp_files: &[PathBuf],
    out_pha: impl AsRef<Path>,
    out_rsp: Option<&Path>,
    add_exposure: bool,
) -> Result<(CztiSpectrum, Option<CztiResponse>)> {
    if pha_files.is_empty() {
        return Err(CztiError::InvalidConfig(
            "run_addspec requires at least one input .pha file".into(),
        ));
    }
    let mut first_meta: Option<CztiObsMeta> = None;
    let mut specs = Vec::with_capacity(pha_files.len());
    for p in pha_files {
        let (m, s) = read_spectrum_file(p)?;
        if first_meta.is_none() {
            first_meta = Some(m);
        }
        specs.push(s);
    }
    let meta = first_meta.unwrap_or_default();

    let rsp_name = out_rsp
        .and_then(|p| p.file_name())
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "NONE".to_string());

    let comb_spec = combine_spectra(&specs, add_exposure, rsp_name);
    write_spectrum_file(&out_pha, &meta, &comb_spec)?;

    let comb_rsp = if let Some(rsp_out) = out_rsp {
        if !rsp_files.is_empty() {
            let mut rsps = Vec::with_capacity(rsp_files.len());
            for rp in rsp_files {
                let (_, r) = read_response_file(rp)?;
                rsps.push(r);
            }
            let exps: Vec<f64> = specs.iter().map(|s| s.exposure).collect();
            let cr = combine_responses(&rsps, &exps, add_exposure);
            write_response_file(rsp_out, &meta, &cr)?;
            Some(cr)
        } else {
            None
        }
    } else {
        None
    };

    Ok((comb_spec, comb_rsp))
}
