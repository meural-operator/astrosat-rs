//! FITS 4.0 zero-copy I/O benchmarks (real AstroSat L1 FITS + synthetic large BINTABLE roundtrip).

use std::fs;
use std::path::Path;

use astrosat_fits::{write_fits, BinTableBuilder, FitsFile, Hdu, HduKind};
use rayon::prelude::*;
use tempfile::tempdir;

use crate::metrics::{measure_in_pool, BenchmarkResult, RunOutcome};

/// Benchmark zero-copy `memmap2` reading & row/column decoding across real AstroSat L1 FITS files.
pub fn bench_fits_real(
    czti_l1_m0: &Path,
    laxpc_l1_ea: &Path,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    if !czti_l1_m0.is_file() {
        anyhow::bail!("missing CZTI L1 file: {}", czti_l1_m0.display());
    }
    if !laxpc_l1_ea.is_file() {
        anyhow::bail!("missing LAXPC L1 file: {}", laxpc_l1_ea.display());
    }

    let czti_bytes = fs::metadata(czti_l1_m0)?.len() as usize;
    let laxpc_bytes = fs::metadata(laxpc_l1_ea)?.len() as usize;

    measure_in_pool(
        "fits_io_real",
        "fits_io",
        "Zero-copy memmap2 open + BINTABLE scan on real CZTI M0 & LAXPC EA L1 FITS",
        threads,
        iterations,
        || {
            let paths = [czti_l1_m0, laxpc_l1_ea];
            let totals: anyhow::Result<Vec<(usize, u64)>> = paths
                .par_iter()
                .map(|path| -> anyhow::Result<(usize, u64)> {
                    let fits = FitsFile::open(path)?;
                    let mut rows = 0usize;
                    let mut acc = 0u64;
                    for hdu in fits.iter_hdus() {
                        if *hdu.kind() == HduKind::BinTable {
                            let bt = hdu.as_bintable()?;
                            rows += bt.nrows();
                            for r in 0..bt.nrows() {
                                let slice = bt.row_bytes(r);
                                if let Some(&first) = slice.first() {
                                    acc = acc.wrapping_add(first as u64);
                                }
                                if slice.len() > 64 {
                                    if let Some(&last) = slice.last() {
                                        acc = acc.wrapping_add(last as u64);
                                    }
                                }
                            }
                        }
                    }
                    Ok((rows, acc))
                })
                .collect();
            let totals = totals?;
            let total_rows: usize = totals.iter().map(|t| t.0).sum();
            let guard_sum: u64 = totals.iter().map(|t| t.1).sum();

            Ok(RunOutcome {
                events_processed: total_rows,
                bytes_processed: czti_bytes + laxpc_bytes,
                notes: format!(
                    "CZTI+LAXPC L1 BINTABLE rows={total_rows}, checksum_guard={guard_sum}"
                ),
            })
        },
    )
}

/// Benchmark synthetic FITS 4.0 `BINTABLE` write (with Seaman/Pence `CHECKSUM`/`DATASUM`)
/// + `memmap2` zero-copy read + typed column decode for `n_rows` event records.
pub fn bench_fits_synth(
    n_rows: usize,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    let dir = tempdir()?;
    let out_path = dir.path().join("synth_events.fits");

    // Pre-build columns outside timing loop so we isolate FITS serialization + I/O + decode.
    let times: Vec<f64> = (0..n_rows).map(|i| 190_000_000.0 + (i as f64) * 1e-4).collect();
    let phas: Vec<u16> = (0..n_rows).map(|i| (i % 1024) as u16).collect();
    let pis: Vec<i16> = (0..n_rows).map(|i| (i % 512) as i16).collect();
    let energies: Vec<f32> = (0..n_rows).map(|i| 5.0 + 0.25 * ((i % 512) as f32)).collect();
    let detx: Vec<u8> = (0..n_rows).map(|i| (i % 64) as u8).collect();
    let dety: Vec<u8> = (0..n_rows).map(|i| ((i / 64) % 64) as u8).collect();

    measure_in_pool(
        "fits_io_synth",
        "fits_io",
        &format!("FITS 4.0 write (CHECKSUM/DATASUM) + memmap2 read ({n_rows} events)"),
        threads,
        iterations,
        || {
            let primary = Hdu::empty_primary();
            let ev_hdu = BinTableBuilder::new("EVENTS")
                .add_f64("TIME", Some("s"), times.clone())
                .add_u16("PHA", Some("chan"), phas.clone())
                .add_i16("PI", Some("chan"), pis.clone())
                .add_f32("ENERGY", Some("keV"), energies.clone())
                .add_u8("DETX", None, detx.clone())
                .add_u8("DETY", None, dety.clone())
                .build()?;

            let mut hdus = vec![primary, ev_hdu];
            write_fits(&out_path, &mut hdus, true)?;
            let file_bytes = fs::metadata(&out_path)?.len() as usize;

            let fits = FitsFile::open(&out_path)?;
            let h1 = fits.hdu_by_name("EVENTS")?;
            let (ds_ok, cs_ok) = h1.verify_checksum();
            anyhow::ensure!(ds_ok && cs_ok, "FITS checksum verification failed");
            let bt = h1.as_bintable()?;
            let t_back = bt.read_col_f64("TIME")?;
            let pi_back = bt.read_col_i16("PI")?;
            let e_back = bt.read_col_f32("ENERGY")?;
            let sum_guard = t_back.len() + pi_back.len() + e_back.len();

            Ok(RunOutcome {
                events_processed: n_rows,
                bytes_processed: file_bytes * 2, // write + read
                notes: format!(
                    "file={:.2} MiB R+W, checksum=OK, cols_decoded={sum_guard}",
                    (file_bytes as f64) / (1024.0 * 1024.0)
                ),
            })
        },
    )
}
