//! Reproducible benchmark suite for `astrosat-rs`.
//!
//! Supports running individual isolated benchmarks (`--suite <target>` via child `--worker`
//! processes) so `/proc/self/status` `VmHWM` peak RSS is uncontaminated across workloads and
//! thread counts, as well as full multi-core scaling sweeps (`1, 4, 16, max` threads)
//! outputting Markdown and JSON.

mod bench_bary;
mod bench_czti;
mod bench_fits;
mod bench_laxpc;
mod metrics;

use std::fs;
use std::path::PathBuf;
use std::process::Command;

use clap::{Parser, ValueEnum};

use metrics::{read_hostname, read_total_ram_gib, BenchmarkReport, BenchmarkResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum SuiteTarget {
    All,
    FitsReal,
    FitsSynth,
    BarySynth,
    BaryFitsReal,
    CztiReal,
    CztiSynth,
    LaxpcReal,
    LaxpcSynth,
}

impl SuiteTarget {
    fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::FitsReal => "fits-real",
            Self::FitsSynth => "fits-synth",
            Self::BarySynth => "bary-synth",
            Self::BaryFitsReal => "bary-fits-real",
            Self::CztiReal => "czti-real",
            Self::CztiSynth => "czti-synth",
            Self::LaxpcReal => "laxpc-real",
            Self::LaxpcSynth => "laxpc-synth",
        }
    }
}

#[derive(Debug, Parser)]
#[command(
    name = "astrosat-bench",
    version,
    about = "Reproducible performance, scaling, throughput, and RSS benchmark suite for astrosat-rs"
)]
struct Cli {
    /// Benchmarks to execute (`all` runs the full suite).
    #[arg(long, value_enum, default_value_t = SuiteTarget::All)]
    suite: SuiteTarget,

    /// Comma-separated thread counts to sweep (e.g. `1,4,16,max`).
    #[arg(long, default_value = "1,4,16,max")]
    threads: String,

    /// Measured iterations per benchmark point.
    #[arg(long, default_value_t = 3)]
    iterations: usize,

    /// Iterations for full real-dataset E2E pipeline runs (czti-real / laxpc-real).
    #[arg(long, default_value_t = 1)]
    real_iterations: usize,

    /// Synthetic FITS / barycentric event count.
    #[arg(long, default_value_t = 1_000_000)]
    synth_events: usize,

    /// Synthetic CZTI total event count across 4 quadrants.
    #[arg(long, default_value_t = 2_000_000)]
    czti_synth_events: usize,

    /// Synthetic LAXPC frame count (2242 bytes/frame).
    #[arg(long, default_value_t = 20_000)]
    laxpc_synth_frames: usize,

    /// Emit machine-readable JSON output to stdout (unless `--json-out` is specified).
    #[arg(long)]
    json: bool,

    /// Write full JSON report to this file path.
    #[arg(long)]
    json_out: Option<PathBuf>,

    /// Write Markdown summary report to this file path.
    #[arg(long)]
    md_out: Option<PathBuf>,

    /// Run in-process instead of spawning isolated per-(benchmark, threads) child processes.
    #[arg(long)]
    in_process: bool,

    /// Internal worker flag: run exactly one target at one thread count and print single JSON row.
    #[arg(long, hide = true)]
    worker: bool,

    /// Sample CZTI Level-1 directory (`9000000406`).
    #[arg(
        long,
        default_value = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1"
    )]
    czti_l1_dir: PathBuf,

    /// Sample CZTI Level-2 reference clean event file (for `bary-fits-real`).
    #[arg(
        long,
        default_value = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0/AS1T01_112T01_9000000406czt_level2_quad_clean.evt"
    )]
    czti_l2_evt: PathBuf,

    /// Sample CZTI orbit file (`9000000406`).
    #[arg(
        long,
        default_value = "/mnt/d/ISRO/data/sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb"
    )]
    czti_orb: PathBuf,

    /// CZTI CALDB directory.
    #[arg(long, default_value = "/mnt/d/ISRO/data/caldb/czti/data/as1/czti")]
    czti_caldb: PathBuf,

    /// Sample LAXPC Level-1 orbit directory (`9000000252` orbit `01509`).
    #[arg(
        long,
        default_value = "/mnt/d/ISRO/data/sample/laxpc_l1/20160107_G02_010T01_9000000252_level1/laxpc/01509"
    )]
    laxpc_l1_dir: PathBuf,

    /// LAXPC CALDB directory.
    #[arg(long, default_value = "/mnt/d/ISRO/data/caldb/laxpc")]
    laxpc_caldb: PathBuf,
}

fn parse_threads(spec: &str, host_cpus: usize) -> anyhow::Result<Vec<usize>> {
    let mut out = Vec::new();
    for tok in spec.split(',') {
        let t = tok.trim();
        if t.is_empty() {
            continue;
        }
        let n = if t.eq_ignore_ascii_case("max") {
            host_cpus.max(1)
        } else {
            let parsed = t
                .parse::<usize>()
                .map_err(|_| anyhow::anyhow!("invalid thread count token: '{t}'"))?;
            if parsed == 0 {
                anyhow::bail!("thread count must be >= 1, got '0'");
            }
            parsed.min(host_cpus.max(1))
        };
        if !out.contains(&n) {
            out.push(n);
        }
    }
    if out.is_empty() {
        out.push(host_cpus.max(1));
    }
    out.sort_unstable();
    Ok(out)
}

fn expand_targets(suite: SuiteTarget) -> Vec<SuiteTarget> {
    match suite {
        SuiteTarget::All => vec![
            SuiteTarget::FitsReal,
            SuiteTarget::FitsSynth,
            SuiteTarget::BarySynth,
            SuiteTarget::BaryFitsReal,
            SuiteTarget::CztiReal,
            SuiteTarget::CztiSynth,
            SuiteTarget::LaxpcReal,
            SuiteTarget::LaxpcSynth,
        ],
        other => vec![other],
    }
}

fn run_single_target(cli: &Cli, target: SuiteTarget, threads: usize) -> anyhow::Result<BenchmarkResult> {
    match target {
        SuiteTarget::FitsReal => {
            let czti_m0 = cli
                .czti_l1_dir
                .join("modeM0/AS1T01_112T01_9000000406cztM0_level1.fits");
            let laxpc_ea = cli
                .laxpc_l1_dir
                .join("lxp1/modeEA/AS1G02_010T01_9000000252lxp1EA_level1.fits");
            bench_fits::bench_fits_real(&czti_m0, &laxpc_ea, threads, cli.iterations)
        }
        SuiteTarget::FitsSynth => {
            bench_fits::bench_fits_synth(cli.synth_events, threads, cli.iterations)
        }
        SuiteTarget::BarySynth => {
            bench_bary::bench_bary_synth(cli.synth_events, threads, cli.iterations)
        }
        SuiteTarget::BaryFitsReal => bench_bary::bench_bary_fits_real(
            &cli.czti_l2_evt,
            &cli.czti_orb,
            83.6331,
            22.0145,
            threads,
            cli.iterations,
        ),
        SuiteTarget::CztiReal => bench_czti::bench_czti_real(
            &cli.czti_l1_dir,
            &cli.czti_caldb,
            threads,
            cli.real_iterations,
        ),
        SuiteTarget::CztiSynth => {
            bench_czti::bench_czti_synth(cli.czti_synth_events, threads, cli.iterations)
        }
        SuiteTarget::LaxpcReal => bench_laxpc::bench_laxpc_real(
            &cli.laxpc_l1_dir,
            Some(&cli.laxpc_caldb),
            threads,
            cli.real_iterations,
        ),
        SuiteTarget::LaxpcSynth => {
            bench_laxpc::bench_laxpc_synth(cli.laxpc_synth_frames, threads, cli.iterations)
        }
        SuiteTarget::All => unreachable!(),
    }
}

fn spawn_worker(cli: &Cli, target: SuiteTarget, threads: usize) -> anyhow::Result<BenchmarkResult> {
    let exe = std::env::current_exe()?;
    let output = Command::new(exe)
        .arg("--worker")
        .arg("--suite")
        .arg(target.as_str())
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--iterations")
        .arg(cli.iterations.to_string())
        .arg("--real-iterations")
        .arg(cli.real_iterations.to_string())
        .arg("--synth-events")
        .arg(cli.synth_events.to_string())
        .arg("--czti-synth-events")
        .arg(cli.czti_synth_events.to_string())
        .arg("--laxpc-synth-frames")
        .arg(cli.laxpc_synth_frames.to_string())
        .arg("--czti-l1-dir")
        .arg(&cli.czti_l1_dir)
        .arg("--czti-l2-evt")
        .arg(&cli.czti_l2_evt)
        .arg("--czti-orb")
        .arg(&cli.czti_orb)
        .arg("--czti-caldb")
        .arg(&cli.czti_caldb)
        .arg("--laxpc-l1-dir")
        .arg(&cli.laxpc_l1_dir)
        .arg("--laxpc-caldb")
        .arg(&cli.laxpc_caldb)
        .output()?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!(
            "worker `{}` ({}T) failed: {}",
            target.as_str(),
            threads,
            stderr.trim()
        );
    }
    let stdout = String::from_utf8(output.stdout)?;
    let row: BenchmarkResult = serde_json::from_str(stdout.trim())?;
    Ok(row)
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    let host_cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let thread_sweep = parse_threads(&cli.threads, host_cpus)?;

    if cli.worker {
        let targets = expand_targets(cli.suite);
        anyhow::ensure!(
            targets.len() == 1 && thread_sweep.len() == 1,
            "--worker expects a single suite target and single thread count"
        );
        let res = run_single_target(&cli, targets[0], thread_sweep[0])?;
        println!("{}", serde_json::to_string(&res)?);
        return Ok(());
    }

    let targets = expand_targets(cli.suite);
    let mut results = Vec::new();
    for target in targets {
        for &threads in &thread_sweep {
            eprintln!(
                "[astrosat-bench] running `{}` with {} thread(s)...",
                target.as_str(),
                threads
            );
            let res = if cli.in_process {
                run_single_target(&cli, target, threads)?
            } else {
                spawn_worker(&cli, target, threads)?
            };
            results.push(res);
        }
    }

    let mut report = BenchmarkReport {
        host_cpus,
        total_ram_gib: read_total_ram_gib(),
        hostname: read_hostname(),
        thread_sweep,
        results,
    };
    report.compute_scaling();

    let md = report.to_markdown();
    let json_str = serde_json::to_string_pretty(&report)?;

    if let Some(ref p) = cli.json_out {
        if let Some(parent) = p.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent)?;
            }
        }
        fs::write(p, &json_str)?;
    }
    if let Some(ref p) = cli.md_out {
        if let Some(parent) = p.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent)?;
            }
        }
        fs::write(p, &md)?;
    }

    if cli.json {
        println!("{json_str}");
    } else {
        println!("{md}");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thread_spec_parses_and_dedups() {
        let parsed = parse_threads("1,4,16,max,4", 52).unwrap();
        assert_eq!(parsed, vec![1, 4, 16, 52]);
        assert!(parse_threads("0", 52).is_err());
        assert!(parse_threads("abc", 52).is_err());
    }

    #[test]
    fn scaling_and_markdown_table_format() {
        let mut report = BenchmarkReport {
            host_cpus: 52,
            total_ram_gib: 62.0,
            hostname: "test-host".into(),
            thread_sweep: vec![1, 4],
            results: vec![
                BenchmarkResult {
                    id: "bary_synth".into(),
                    category: "bary".into(),
                    description: "test".into(),
                    threads: 1,
                    iterations: 3,
                    events_processed: 100_000,
                    bytes_processed: 1_600_000,
                    wall_ms: 100.0,
                    latency: metrics::LatencyStats::from_samples_ms(&[99.0, 100.0, 101.0]),
                    events_per_sec: 1_000_000.0,
                    mb_per_sec: 16.0,
                    peak_rss_mib: 24.5,
                    speedup_vs_1t: None,
                    notes: "ok".into(),
                },
                BenchmarkResult {
                    id: "bary_synth".into(),
                    category: "bary".into(),
                    description: "test".into(),
                    threads: 4,
                    iterations: 3,
                    events_processed: 100_000,
                    bytes_processed: 1_600_000,
                    wall_ms: 28.0,
                    latency: metrics::LatencyStats::from_samples_ms(&[27.0, 28.0, 29.0]),
                    events_per_sec: 3_571_428.0,
                    mb_per_sec: 57.1,
                    peak_rss_mib: 28.0,
                    speedup_vs_1t: None,
                    notes: "ok".into(),
                },
            ],
        };
        report.compute_scaling();
        assert!((report.results[1].speedup_vs_1t.unwrap() - (100.0 / 28.0)).abs() < 1e-6);
        let md = report.to_markdown();
        assert!(md.contains("bary_synth"));
        assert!(md.contains("3.57x"));
    }

    #[test]
    fn fits_and_bary_and_synth_pipelines_smoke() {
        let r_fits = bench_fits::bench_fits_synth(10_000, 2, 1).unwrap();
        assert_eq!(r_fits.events_processed, 10_000);
        assert!(r_fits.wall_ms > 0.0);

        let r_bary = bench_bary::bench_bary_synth(10_000, 2, 1).unwrap();
        assert_eq!(r_bary.events_processed, 10_000);

        let r_czti = bench_czti::bench_czti_synth(16_000, 2, 1).unwrap();
        assert_eq!(r_czti.events_processed, 16_000);

        let r_laxpc = bench_laxpc::bench_laxpc_synth(1024, 2, 1).unwrap();
        assert!(r_laxpc.events_processed > 0);
    }
}
