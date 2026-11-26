//! Machine-readable (`--json`) and human-readable output formatting utilities.

use serde::Serialize;

/// Global output mode shared across all CLI subcommands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OutputMode {
    pub json: bool,
    pub quiet: bool,
}

impl OutputMode {
    pub fn emit<T: Serialize>(&self, payload: &T, human_fn: impl FnOnce(&T)) -> anyhow::Result<()> {
        if self.json {
            println!("{}", serde_json::to_string_pretty(payload)?);
        } else if !self.quiet {
            human_fn(payload);
        }
        Ok(())
    }
}

/// Render an actionable error message to stderr (or JSON to stdout when `--json` is active).
pub fn format_error_hint(err: &anyhow::Error) -> Option<&'static str> {
    let msg = format!("{err:#}").to_ascii_lowercase();
    if msg.contains("caldb") || msg.contains("caldb.indx") {
        Some("Hint: Pass --caldb <path> pointing to the CZTI or LAXPC CALDB directory (e.g. <caldb_dir>/czti or <caldb_dir>/laxpc).")
    } else if msg.contains("ra_obj") || msg.contains("dec_obj") || msg.contains("ra/dec") {
        Some("Hint: Provide target coordinates explicitly with --ra <deg> --dec <deg>.")
    } else if msg.contains("orbit") || msg.contains(".orb") {
        Some("Hint: Supply an AstroSat .orb FITS file with --orbit <path> (or --mkf <path> to synthesize an orbit from MKF).")
    } else if msg.contains("no orbits") || msg.contains("missing input") || msg.contains("discover") {
        Some("Hint: Pass --input <l1_dir> pointing to an AstroSat Level-1 observation or orbit directory, or specify individual input files explicitly.")
    } else if msg.contains("fits") && (msg.contains("no such file") || msg.contains("cannot find")) {
        Some("Hint: Verify that the input FITS file path exists and is readable.")
    } else {
        None
    }
}

/// Render a compact ASCII bar sparkline for a slice of non-negative values.
pub fn ascii_sparkline(values: &[f64]) -> String {
    const BARS: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
    if values.is_empty() {
        return String::new();
    }
    let min_v = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max_v = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let span = (max_v - min_v).max(1e-12);
    values
        .iter()
        .map(|&v| {
            let idx = (((v - min_v) / span) * 7.0).round().clamp(0.0, 7.0) as usize;
            BARS[idx]
        })
        .collect()
}
