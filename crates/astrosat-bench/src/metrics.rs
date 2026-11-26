//! Timing, throughput, multi-core scaling, and peak RSS memory measurement helpers.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::time::Instant;

/// Memory stats parsed from `/proc/self/status` on Linux (or zeros on non-Linux).
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct ProcMemStats {
    /// Peak resident set size (`VmHWM`) in MiB (`1024 * 1024` bytes).
    pub peak_rss_mib: f64,
    /// Current resident set size (`VmRSS`) in MiB.
    pub current_rss_mib: f64,
}

/// Read `/proc/self/status` to extract `VmHWM` and `VmRSS` in MiB.
pub fn read_proc_mem_stats() -> ProcMemStats {
    let Ok(status) = fs::read_to_string("/proc/self/status") else {
        return ProcMemStats::default();
    };
    let mut stats = ProcMemStats::default();
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            if let Some(kb) = parse_proc_kb(rest) {
                stats.peak_rss_mib = kb / 1024.0;
            }
        } else if let Some(rest) = line.strip_prefix("VmRSS:") {
            if let Some(kb) = parse_proc_kb(rest) {
                stats.current_rss_mib = kb / 1024.0;
            }
        }
    }
    stats
}

fn parse_proc_kb(s: &str) -> Option<f64> {
    s.split_whitespace().next()?.parse::<f64>().ok()
}

/// Read total physical system RAM in GiB from `/proc/meminfo` (or `0.0` if unavailable).
pub fn read_total_ram_gib() -> f64 {
    let Ok(info) = fs::read_to_string("/proc/meminfo") else {
        return 0.0;
    };
    for line in info.lines() {
        if let Some(rest) = line.strip_prefix("MemTotal:") {
            if let Some(kb) = parse_proc_kb(rest) {
                return kb / (1024.0 * 1024.0);
            }
        }
    }
    0.0
}

/// Read system hostname from `/etc/hostname` or fallback.
pub fn read_hostname() -> String {
    fs::read_to_string("/etc/hostname")
        .ok()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "localhost".to_string())
}

/// Latency summary computed across one or more measured iterations.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct LatencyStats {
    pub mean_ms: f64,
    pub min_ms: f64,
    pub median_ms: f64,
    pub p95_ms: f64,
    pub max_ms: f64,
}

impl LatencyStats {
    pub fn from_samples_ms(samples: &[f64]) -> Self {
        if samples.is_empty() {
            return Self {
                mean_ms: 0.0,
                min_ms: 0.0,
                median_ms: 0.0,
                p95_ms: 0.0,
                max_ms: 0.0,
            };
        }
        let mut sorted = samples.to_vec();
        sorted.sort_by(|a, b| a.total_cmp(b));
        let n = sorted.len();
        let mean_ms = sorted.iter().sum::<f64>() / (n as f64);
        let min_ms = sorted[0];
        let max_ms = sorted[n - 1];
        let median_ms = if n % 2 == 1 {
            sorted[n / 2]
        } else {
            0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
        };
        let p95_idx = (((n as f64) * 0.95).ceil() as usize)
            .saturating_sub(1)
            .min(n - 1);
        let p95_ms = sorted[p95_idx];
        Self {
            mean_ms,
            min_ms,
            median_ms,
            p95_ms,
            max_ms,
        }
    }
}

/// Single measured benchmark record for one workload and thread count.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    pub id: String,
    pub category: String,
    pub description: String,
    pub threads: usize,
    pub iterations: usize,
    pub events_processed: usize,
    pub bytes_processed: usize,
    pub wall_ms: f64,
    pub latency: LatencyStats,
    pub events_per_sec: f64,
    pub mb_per_sec: f64,
    pub peak_rss_mib: f64,
    pub speedup_vs_1t: Option<f64>,
    pub notes: String,
}

/// Output produced by one single-iteration/multi-iteration workload closure.
#[derive(Debug, Clone)]
pub struct RunOutcome {
    pub events_processed: usize,
    pub bytes_processed: usize,
    pub notes: String,
}

/// Time a workload `iterations` times inside a dedicated `rayon` pool with `threads` threads.
pub fn measure_in_pool<F>(
    id: &str,
    category: &str,
    description: &str,
    threads: usize,
    iterations: usize,
    mut work: F,
) -> anyhow::Result<BenchmarkResult>
where
    F: FnMut() -> anyhow::Result<RunOutcome> + Send,
{
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()?;

    let mut samples_ms = Vec::with_capacity(iterations.max(1));
    let mut last_outcome = RunOutcome {
        events_processed: 0,
        bytes_processed: 0,
        notes: String::new(),
    };

    for _ in 0..iterations.max(1) {
        let t0 = Instant::now();
        last_outcome = pool.install(&mut work)?;
        let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
        samples_ms.push(elapsed_ms);
    }

    let latency = LatencyStats::from_samples_ms(&samples_ms);
    let wall_ms = latency.mean_ms;
    let sec = (wall_ms / 1000.0).max(1e-9);
    let events_per_sec = (last_outcome.events_processed as f64) / sec;
    let mb_per_sec = (last_outcome.bytes_processed as f64) / (1_000_000.0 * sec);
    let mem = read_proc_mem_stats();

    Ok(BenchmarkResult {
        id: id.to_string(),
        category: category.to_string(),
        description: description.to_string(),
        threads: threads.max(1),
        iterations: iterations.max(1),
        events_processed: last_outcome.events_processed,
        bytes_processed: last_outcome.bytes_processed,
        wall_ms,
        latency,
        events_per_sec,
        mb_per_sec,
        peak_rss_mib: mem.peak_rss_mib,
        speedup_vs_1t: None,
        notes: last_outcome.notes,
    })
}

/// Complete benchmark suite report with system metadata and scaling calculations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkReport {
    pub host_cpus: usize,
    pub total_ram_gib: f64,
    pub hostname: String,
    pub thread_sweep: Vec<usize>,
    pub results: Vec<BenchmarkResult>,
}

impl BenchmarkReport {
    pub fn compute_scaling(&mut self) {
        let mut baselines: BTreeMap<String, f64> = BTreeMap::new();
        for r in &self.results {
            if r.threads == 1 && r.wall_ms > 0.0 {
                baselines.insert(r.id.clone(), r.wall_ms);
            }
        }
        for r in &mut self.results {
            if let Some(&base_ms) = baselines.get(&r.id) {
                r.speedup_vs_1t = Some(base_ms / r.wall_ms.max(1e-9));
            }
        }
    }

    pub fn to_markdown(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "# `astrosat-rs` Benchmark Report\n\n- **Host**: `{}` ({} vCPUs, {:.1} GiB RAM)\n- **Thread Sweep**: `{:?}`\n\n",
            self.hostname, self.host_cpus, self.total_ram_gib, self.thread_sweep
        ));

        out.push_str("## Summary Results\n\n");
        out.push_str("| Benchmark | Threads | Wall (ms) | Median (ms) | Events/s | I/O (MB/s) | Speedup vs 1T | Peak RSS (MiB) | Notes |\n");
        out.push_str("|---|---:|---:|---:|---:|---:|---:|---:|---|\n");
        for r in &self.results {
            let speedup = r
                .speedup_vs_1t
                .map(|s| format!("{s:.2}x"))
                .unwrap_or_else(|| "-".to_string());
            let ev_s = if r.events_processed > 0 {
                format!("{:.0}", r.events_per_sec)
            } else {
                "-".to_string()
            };
            let mb_s = if r.bytes_processed > 0 {
                format!("{:.1}", r.mb_per_sec)
            } else {
                "-".to_string()
            };
            out.push_str(&format!(
                "| `{}` | {} | {:.2} | {:.2} | {} | {} | {} | {:.1} | {} |\n",
                r.id,
                r.threads,
                r.wall_ms,
                r.latency.median_ms,
                ev_s,
                mb_s,
                speedup,
                r.peak_rss_mib,
                r.notes
            ));
        }

        // Group scaling summary by benchmark ID
        let mut by_id: BTreeMap<&str, Vec<&BenchmarkResult>> = BTreeMap::new();
        for r in &self.results {
            by_id.entry(r.id.as_str()).or_default().push(r);
        }
        out.push_str("\n## Multi-Core Scaling Summary\n\n");
        out.push_str("| Benchmark | 1T (ms) | 4T (ms) | 16T (ms) | MaxT (ms) | Best Speedup | Peak RSS @ MaxT (MiB) |\n");
        out.push_str("|---|---:|---:|---:|---:|---:|---:|\n");
        for (id, mut rows) in by_id {
            rows.sort_by_key(|r| r.threads);
            let find_t = |t: usize| -> Option<&BenchmarkResult> {
                rows.iter().copied().find(|r| r.threads == t)
            };
            let max_row = rows.last().copied();
            let fmt_ms = |opt: Option<&BenchmarkResult>| -> String {
                opt.map(|r| format!("{:.2}", r.wall_ms))
                    .unwrap_or_else(|| "-".to_string())
            };
            let best_speedup = rows
                .iter()
                .filter_map(|r| r.speedup_vs_1t)
                .fold(None::<f64>, |acc, v| Some(acc.map_or(v, |a| a.max(v))))
                .map(|s| format!("{s:.2}x"))
                .unwrap_or_else(|| "-".to_string());
            let max_rss = max_row
                .map(|r| format!("{:.1}", r.peak_rss_mib))
                .unwrap_or_else(|| "-".to_string());
            out.push_str(&format!(
                "| `{}` | {} | {} | {} | {} | {} | {} |\n",
                id,
                fmt_ms(find_t(1)),
                fmt_ms(find_t(4)),
                fmt_ms(find_t(16)),
                fmt_ms(max_row),
                best_speedup,
                max_rss
            ));
        }
        out
    }
}
