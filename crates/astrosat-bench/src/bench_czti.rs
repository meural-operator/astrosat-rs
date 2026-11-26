//! AstroSat CZTI Level-1 -> Level-2 benchmarks (real ObsID 9000000406 + synthetic large stream).

use std::fs;
use std::path::Path;

use astrosat_czti::fits_io::{write_l2_event_file, EventFileWriteOpts};
use astrosat_czti::{
    bunchclean_quad, evtclean_quad, evtencor_quad, pha2energy_quad, run_pipeline, BunchCleanConfig,
    BunchTable, CztiGtiOutput, CztiObsMeta, CztiPipelineConfig, EboundsTable, EvtCleanConfig,
    EvtEnCorConfig, QuadEvents, QuadGainCorCaldb, QuadGainOffset, QuadTempSeries, SsmDataTable,
    TempTable, VetoSpectrumTable, NUM_QUADS, PIXELS_PER_QUAD,
};
use rayon::prelude::*;
use tempfile::tempdir;

use crate::metrics::{measure_in_pool, BenchmarkResult, RunOutcome};

fn dir_files_bytes(dir: &Path) -> usize {
    let Ok(entries) = fs::read_dir(dir) else {
        return 0;
    };
    entries
        .filter_map(|e| e.ok())
        .filter_map(|e| e.metadata().ok())
        .filter(|m| m.is_file())
        .map(|m| m.len() as usize)
        .sum()
}

/// Benchmark full end-to-end CZTI L1 -> L2 pipeline (`run_pipeline`) on real ObsID `9000000406`.
pub fn bench_czti_real(
    l1_dir: &Path,
    caldb_dir: &Path,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    let l1_m0 = l1_dir.join("modeM0/AS1T01_112T01_9000000406cztM0_level1.fits");
    let l1_ss = l1_dir.join("modeSS/AS1T01_112T01_9000000406cztSS_level1.fits");
    let tct = l1_dir.join("aux/AS1T01_112T01_9000000406czt_level1.tct");
    let mkf = l1_dir.join("AS1T01_112T01_9000000406czt_level2.mkf");
    let mkf_thresh = l1_dir.join("mkfThresholds.txt");

    for p in [&l1_m0, &tct, &mkf] {
        if !p.is_file() {
            anyhow::bail!("missing CZTI L1 input: {}", p.display());
        }
    }
    if !caldb_dir.exists() {
        anyhow::bail!("missing CZTI CALDB dir: {}", caldb_dir.display());
    }

    let in_bytes = fs::metadata(&l1_m0)?.len() as usize
        + fs::metadata(&l1_ss).map(|m| m.len() as usize).unwrap_or(0)
        + fs::metadata(&tct)?.len() as usize
        + fs::metadata(&mkf)?.len() as usize;

    measure_in_pool(
        "czti_real_9000000406",
        "czti_real",
        "Full CZTI L1->L2 pipeline (Stages 1-3) on real ObsID 9000000406",
        threads,
        iterations,
        || {
            let out_dir = tempdir()?;
            let mut cfg = CztiPipelineConfig::new(
                &l1_m0,
                l1_ss.is_file().then_some(l1_ss.clone()),
                &tct,
                &mkf,
                caldb_dir,
                out_dir.path(),
            );
            if mkf_thresh.is_file() {
                cfg.mkf_threshold_file = Some(mkf_thresh.clone());
            }
            let summary = run_pipeline(&cfg)?;
            let raw_total: usize = summary.stage1_raw_events.iter().sum();
            let clean_total: usize = summary.stage2_clean_events.iter().sum();
            let out_bytes = dir_files_bytes(out_dir.path());
            Ok(RunOutcome {
                events_processed: raw_total,
                bytes_processed: in_bytes + out_bytes,
                notes: format!(
                    "ObsID 9000000406: raw={raw_total}, clean={clean_total}, exp_q0={:.1}s",
                    summary.quad_exposures_sec[0]
                ),
            })
        },
    )
}

/// Build `n_events` realistic synthetic CZTI events for one segment of quadrant `quad_id`.
fn make_synth_quad_segment(quad_id: u8, n_events: usize, t_start: f64) -> QuadEvents {
    let mut ev = QuadEvents::with_capacity(quad_id, n_events);
    ev.tstart = t_start;
    ev.tstop = t_start + (n_events as f64) * 5e-4;
    ev.exposure = ev.tstop - ev.tstart;
    for i in 0..n_events {
        let t = t_start + (i as f64) * 5e-4;
        let det = ((i * 7) % 16) as u8;
        let pix = ((i * 13) % 256) as u8;
        let (dx, dy) = astrosat_czti::detid_pixid_to_detxy(quad_id, det, pix);
        ev.time.push(t);
        ev.cztseccnt.push(t.floor());
        ev.cztntick.push(((t.fract()) * 50_000.0) as u16);
        ev.pha.push(40 + ((i * 19) % 450) as u16);
        ev.detid.push(det);
        ev.pixid.push(pix);
        ev.detx.push(dx);
        ev.dety.push(dy);
        ev.veto.push(if i % 97 == 0 { 10 } else { 0 });
        ev.alpha.push(if i % 131 == 0 { 1 } else { 0 });
    }
    ev
}

/// Benchmark synthetic large CZTI event stream (`n_events_total` across 16 parallel segments & 4 quadrants)
/// through `bunchclean -> pha2energy -> evtencor -> evtclean -> FITS write`.
pub fn bench_czti_synth(
    n_events_total: usize,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    const N_SEGS: usize = 16;
    let per_seg = (n_events_total / N_SEGS).max(1000);
    let total_events = per_seg * N_SEGS;
    let seg_dur = (per_seg as f64) * 5e-4 + 1.0;

    let segments: Vec<(usize, QuadEvents)> = (0..N_SEGS)
        .map(|s| {
            let q = s % NUM_QUADS;
            let t0 = 195_000_000.0 + (s as f64) * seg_dur;
            (q, make_synth_quad_segment(q as u8, per_seg, t0))
        })
        .collect();

    let ebounds = EboundsTable::default();
    let gain_off = QuadGainOffset::default();
    let gain_caldb = QuadGainCorCaldb::default();
    let temp_series = QuadTempSeries {
        time: vec![195_000_000.0, 195_050_000.0],
        temp: vec![[10.0; 16], [11.0; 16]],
    };
    let bc_cfg = BunchCleanConfig::default();
    let encor_cfg = EvtEnCorConfig::default();
    let clean_cfg = EvtCleanConfig::default();

    measure_in_pool(
        "czti_synth",
        "czti_synth",
        &format!("Synthetic large CZTI stream ({total_events} events, 16 segs: bunchclean+pha2energy+evtencor+evtclean+FITS)"),
        threads,
        iterations,
        || {
            let processed_segs: Vec<(usize, QuadEvents)> = segments
                .par_iter()
                .map(|(q, ev)| {
                    let gti = vec![(ev.tstart, ev.tstop)];
                    let bunches = BunchTable::new(*q as u8);
                    let mut bc = bunchclean_quad(ev, &bunches, &gti, &bc_cfg).events;
                    pha2energy_quad(&mut bc, &temp_series, &gain_off, &ebounds);
                    let ec = evtencor_quad(&bc, &gain_caldb, ev.tstart, false, &encor_cfg);
                    let cleaned = evtclean_quad(&ec.events, &clean_cfg);
                    (*q, cleaned)
                })
                .collect();

            // Merge back into 4 quadrant event tables and write L2 FITS
            let mut quads = [
                QuadEvents::new(0),
                QuadEvents::new(1),
                QuadEvents::new(2),
                QuadEvents::new(3),
            ];
            let mut clean_total = 0usize;
            for (q, seg) in processed_segs {
                clean_total += seg.len();
                let dst = &mut quads[q];
                dst.time.extend_from_slice(&seg.time);
                dst.cztseccnt.extend_from_slice(&seg.cztseccnt);
                dst.cztntick.extend_from_slice(&seg.cztntick);
                dst.pha.extend_from_slice(&seg.pha);
                dst.detid.extend_from_slice(&seg.detid);
                dst.pixid.extend_from_slice(&seg.pixid);
                dst.detx.extend_from_slice(&seg.detx);
                dst.dety.extend_from_slice(&seg.dety);
                dst.veto.extend_from_slice(&seg.veto);
                dst.alpha.extend_from_slice(&seg.alpha);
                dst.pi.extend_from_slice(&seg.pi);
                dst.energy.extend_from_slice(&seg.energy);
                dst.exposure += seg.exposure;
            }

            let out_dir = tempdir()?;
            let out_evt = out_dir.path().join("synth_czti_clean.evt");
            let meta = CztiObsMeta::default();
            let gtis = CztiGtiOutput::default();
            let opts = EventFileWriteOpts {
                pixel_exposure: Some(std::array::from_fn(|_| vec![1.0; PIXELS_PER_QUAD])),
                ..Default::default()
            };
            write_l2_event_file(
                &out_evt,
                &meta,
                &quads,
                &VetoSpectrumTable::default(),
                &SsmDataTable::default(),
                &TempTable::default(),
                &gtis,
                &opts,
            )?;
            let out_bytes = fs::metadata(&out_evt)?.len() as usize;

            Ok(RunOutcome {
                events_processed: total_events,
                bytes_processed: out_bytes,
                notes: format!(
                    "raw={total_events}, clean={clean_total}, evt_fits={:.2} MiB",
                    (out_bytes as f64) / (1024.0 * 1024.0)
                ),
            })
        },
    )
}
