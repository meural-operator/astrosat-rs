//! AstroSat LAXPC Level-1 -> Level-2 benchmarks (real Crab ObsID 9000000252 + synthetic large frame stream).

use std::fs;
use std::path::Path;

use astrosat_fits::OgipMeta;
use astrosat_laxpc::{
    bin_events_raw, build_calibrated_events, decode_ea_orbits, run_pipeline, sort_and_dedup_frames,
    write_format_b_event_fits, EventFilter, FrameMode, LaxpcUnit, PipelineConfig,
    COARSE_TIME_MARKER, END_OF_FRAME_MARKER, FRAME_SYNC_BYTE, L1_FRAME_SIZE, MODE_EA_BYTE,
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

/// Benchmark full multi-unit LAXPC L1 -> L2 pipeline (`run_pipeline`) on Crab ObsID `9000000252`.
pub fn bench_laxpc_real(
    l1_dir: &Path,
    caldb_dir: Option<&Path>,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    if !l1_dir.exists() {
        anyhow::bail!("missing LAXPC L1 dir: {}", l1_dir.display());
    }

    // Sum L1 FITS/MKF/TCT sizes under l1_dir
    let mut in_bytes = 0usize;
    for unit in LaxpcUnit::ALL {
        let p = unit.prefix();
        in_bytes += dir_files_bytes(&l1_dir.join(p).join("modeEA"));
        in_bytes += dir_files_bytes(&l1_dir.join(p).join("modeBB"));
    }
    in_bytes += dir_files_bytes(&l1_dir.join("aux"));
    in_bytes += dir_files_bytes(l1_dir);

    measure_in_pool(
        "laxpc_real_9000000252",
        "laxpc_real",
        "Full LAXPC10/20/30 L1->L2 pipeline (Format B + Format A + Orbit) on Crab ObsID 9000000252",
        threads,
        iterations,
        || {
            let out_dir = tempdir()?;
            let cfg = PipelineConfig {
                input_dir: Some(l1_dir.to_path_buf()),
                output_dir: out_dir.path().to_path_buf(),
                caldb_dir: caldb_dir.filter(|p| p.exists()).map(Path::to_path_buf),
                write_format_b_events: true,
                write_format_a_combined: true,
                write_orbit_file: true,
                ..Default::default()
            };
            let res = run_pipeline(&cfg)?;
            let valid_total: usize = res.unit_results.iter().map(|u| u.total_valid_events).sum();
            let gti_total: usize = res.unit_results.iter().map(|u| u.gti_filtered_events).sum();
            let out_bytes = dir_files_bytes(out_dir.path());
            Ok(RunOutcome {
                events_processed: valid_total,
                bytes_processed: in_bytes + out_bytes,
                notes: format!(
                    "ObsID 9000000252 (Lxp10/20/30): valid={valid_total}, gti_filt={gti_total}, fmtA={}",
                    res.format_a_total_events
                ),
            })
        },
    )
}

/// Build `n_frames` valid 2242-byte Level-1 EA telemetry frames with 200 valid events per frame.
fn build_synthetic_ea_frames(n_frames: usize, unit: LaxpcUnit) -> Vec<u8> {
    let mut buf = vec![0u8; n_frames * L1_FRAME_SIZE];
    // Align to 256-tick coarse boundary so each frame has stbg = base + f * 256 (m2 advances by +1)
    let base_stbg: u64 = 1_000_000_000 & !0xFF;

    for f in 0..n_frames {
        let frame = &mut buf[f * L1_FRAME_SIZE..(f + 1) * L1_FRAME_SIZE];
        let stbg = base_stbg + (f as u64) * 256;
        let l1_time = 181_218_232.191_64 + (stbg as f64) / 100_000.0;
        frame[0..8].copy_from_slice(&l1_time.to_be_bytes());
        frame[36..44].copy_from_slice(&(f as i64).to_be_bytes());
        frame[44] = FRAME_SYNC_BYTE;
        frame[45] = unit.id();
        frame[46] = MODE_EA_BYTE;

        let stbg_be = stbg.to_be_bytes();
        frame[49..54].copy_from_slice(&stbg_be[3..8]);
        frame[57] = 0; // valid status_flag
        let ifr = (f % 65536) as u16;
        frame[58..60].copy_from_slice(&ifr.to_be_bytes());

        // Write 1 coarse-time marker (0xEF, m4, m3, m2, 0) followed by 200 single-anode events
        let mut pos = 60usize;
        frame[pos] = COARSE_TIME_MARKER;
        frame[pos + 1] = stbg_be[4]; // m4
        frame[pos + 2] = stbg_be[5]; // m3
        frame[pos + 3] = stbg_be[6]; // m2
        frame[pos + 4] = 0;
        pos += 5;

        for e in 0..200usize {
            if pos + 5 > 2092 {
                break;
            }
            let fine = ((e + 2) & 0xFF) as u8;
            let ia1 = ((e % 7) + 1) as u8; // anodes 1..=7
            let pha: u16 = 50 + ((e * 3 + f) % 400) as u16;
            let ic1_hi = (pha >> 2) as u8;
            let ic1_lo = ((pha & 0x03) << 2) as u8;

            frame[pos] = ia1; // ia2=0, ia1=1..7
            frame[pos + 1] = fine;
            frame[pos + 2] = ic1_hi;
            frame[pos + 3] = 0; // ic2=0 -> single event
            frame[pos + 4] = ic1_lo;
            pos += 5;
        }
        if pos + 5 <= 2092 {
            frame[pos] = END_OF_FRAME_MARKER;
        }
    }
    buf
}

/// Benchmark synthetic large LAXPC EA telemetry stream (`n_frames` x 2242-byte frames across 16 orbit shards):
/// `sort_and_dedup_frames -> decode_ea_orbits -> bin_events_raw -> build_calibrated_events -> Format B FITS write`.
pub fn bench_laxpc_synth(
    n_frames: usize,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    const N_SHARDS: usize = 16;
    let frames_per_shard = (n_frames / N_SHARDS).max(64);
    let total_frames = frames_per_shard * N_SHARDS;
    let raw_bytes = build_synthetic_ea_frames(total_frames, LaxpcUnit::Lxp10);
    let shard_bytes = frames_per_shard * L1_FRAME_SIZE;
    let raw_len = raw_bytes.len();

    measure_in_pool(
        "laxpc_synth",
        "laxpc_synth",
        &format!(
            "Synthetic large LAXPC stream ({total_frames} frames, {:.1} MiB raw: sort+decode+bin+cal+FITS)",
            (raw_len as f64) / (1024.0 * 1024.0)
        ),
        threads,
        iterations,
        || {
            let ordered_orbits: Vec<_> = (0..N_SHARDS)
                .into_par_iter()
                .map(|s| {
                    let slice = &raw_bytes[s * shard_bytes..(s + 1) * shard_bytes];
                    let (chunks, _) = slice.as_chunks::<L1_FRAME_SIZE>();
                    let iter = chunks.iter().map(|c| c.as_slice());
                    sort_and_dedup_frames(iter, Some(LaxpcUnit::Lxp10), Some(FrameMode::Ea))
                        .map(|(frames, _diag)| frames)
                })
                .collect::<astrosat_laxpc::Result<Vec<_>>>()?;

            let stream = decode_ea_orbits(&ordered_orbits, LaxpcUnit::Lxp10);
            let filter = EventFilter::default();
            let _raw_bins = bin_events_raw(&stream.events, 1.0, &filter)?;
            let gti_ins = vec![(stream.t_ins_min - 1.0, stream.t_ins_max + 1.0)];
            let cal = build_calibrated_events(
                &stream.events,
                &gti_ins,
                &filter,
                512,
                |t_ins| t_ins + 190_000_000.0,
            );

            let out_dir = tempdir()?;
            let evn_path = out_dir.path().join("synth_lxp1level2.evn");
            let meta = OgipMeta {
                instrume: "LAXPC".to_string(),
                detnam: Some("LAXPC10".to_string()),
                ..Default::default()
            };
            let gti_utc = vec![(
                stream.t_ins_min + 190_000_000.0,
                stream.t_ins_max + 190_000_000.0,
            )];
            write_format_b_event_fits(
                &evn_path,
                &cal,
                &gti_utc,
                LaxpcUnit::Lxp10,
                &filter,
                512,
                &meta,
            )?;
            let out_bytes = fs::metadata(&evn_path)?.len() as usize;

            Ok(RunOutcome {
                events_processed: stream.total_valid_events,
                bytes_processed: raw_len + out_bytes,
                notes: format!(
                    "frames={total_frames}, decoded_ev={}, out_evn={:.2} MiB",
                    stream.total_valid_events,
                    (out_bytes as f64) / (1024.0 * 1024.0)
                ),
            })
        },
    )
}
