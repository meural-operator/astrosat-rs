//! Solar System Barycentric correction throughput benchmarks.

use std::f64::consts::PI;
use std::fs;
use std::path::Path;

use astrosat_core::{
    barycenter_fits, mjd_utc_to_astrosat_utc, BaryCalculator, BaryEphemeris, BaryFitsOptions,
    CoordFrame, EphemerisKind, OrbitTable,
};
use tempfile::tempdir;

use crate::metrics::{measure_in_pool, BenchmarkResult, RunOutcome};

fn synthetic_leo_orbit() -> anyhow::Result<OrbitTable> {
    let r_km = 7028.0f64;
    let v_kms = 7.53f64;
    let omega = v_kms / r_km;
    let t0 = mjd_utc_to_astrosat_utc(57394.0) - 100.0;
    let n = 1000usize;
    let mut times = Vec::with_capacity(n);
    let mut pos_km = Vec::with_capacity(n);
    let mut vel_kms = Vec::with_capacity(n);
    for i in 0..n {
        let dt = (i as f64) * 6.0;
        let phi = omega * dt + 0.1 * PI;
        times.push(t0 + dt);
        pos_km.push([r_km * phi.cos(), r_km * phi.sin(), 120.0 * (2.0 * phi).sin()]);
        vel_kms.push([
            -v_kms * phi.sin(),
            v_kms * phi.cos(),
            240.0 * omega * (2.0 * phi).cos(),
        ]);
    }
    Ok(OrbitTable::from_arrays(times, pos_km, vel_kms)?)
}

/// Benchmark parallel `BaryCalculator::correct_times` throughput across `n_events`.
pub fn bench_bary_synth(
    n_events: usize,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    let orbit = synthetic_leo_orbit()?;
    let calc = BaryCalculator::new(
        83.6331,
        22.0145,
        CoordFrame::Icrs,
        orbit,
        BaryEphemeris::Analytical,
    );
    let t0 = mjd_utc_to_astrosat_utc(57394.0);
    let times: Vec<f64> = (0..n_events).map(|i| t0 + (i as f64) * 0.004).collect();

    measure_in_pool(
        "bary_synth",
        "bary",
        &format!("Parallel relativistic SSB barycentric correction ({n_events} timestamps)"),
        threads,
        iterations,
        || {
            let tdb = calc.correct_times(&times)?;
            let first = tdb.first().copied().unwrap_or(0.0);
            let last = tdb.last().copied().unwrap_or(0.0);
            Ok(RunOutcome {
                events_processed: tdb.len(),
                bytes_processed: tdb.len() * std::mem::size_of::<f64>() * 2,
                notes: format!(
                    "analytical ephem, TDB span=[{first:.3}, {last:.3}]"
                ),
            })
        },
    )
}

/// Benchmark end-to-end `barycenter_fits` on a real AstroSat CZTI Level-2 event file + `.orb`.
pub fn bench_bary_fits_real(
    evt_fits: &Path,
    orb_fits: &Path,
    ra_deg: f64,
    dec_deg: f64,
    threads: usize,
    iterations: usize,
) -> anyhow::Result<BenchmarkResult> {
    if !evt_fits.is_file() {
        anyhow::bail!("missing event FITS for bary bench: {}", evt_fits.display());
    }
    if !orb_fits.is_file() {
        anyhow::bail!("missing orbit FITS for bary bench: {}", orb_fits.display());
    }

    let in_bytes = fs::metadata(evt_fits)?.len() as usize + fs::metadata(orb_fits)?.len() as usize;
    let dir = tempdir()?;
    let out_path = dir.path().join("bary_out.evt");
    let opts = BaryFitsOptions {
        frame: CoordFrame::Icrs,
        ephemeris: EphemerisKind::Analytical,
        write_checksums: true,
    };

    measure_in_pool(
        "bary_fits_real",
        "bary",
        "End-to-end barycenter_fits on real ObsID 9000000406 CZTI clean .evt + .orb",
        threads,
        iterations,
        || {
            let summary =
                barycenter_fits(evt_fits, &out_path, orb_fits, ra_deg, dec_deg, &opts)?;
            let out_bytes = fs::metadata(&out_path)?.len() as usize;
            Ok(RunOutcome {
                events_processed: summary.events_corrected,
                bytes_processed: in_bytes + out_bytes,
                notes: format!(
                    "hdus={}, gti_rows={}, plephem={}",
                    summary.hdus_updated, summary.gti_rows_corrected, summary.plephem
                ),
            })
        },
    )
}
