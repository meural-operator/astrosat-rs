//! Shared configuration schema (`AppConfig`) and lightweight pure-Rust TOML subset parser
//! (`tables [section]`, nested dotted keys, `#` comments, strings, integers, floats, booleans,
//! and inline arrays) used by both `astrosat-cli --config` and `astrosat-ui` (`POST /api/config/parse`).

use std::path::PathBuf;

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

use crate::error::{CoreError, Result};

/// Top-level CLI configuration loaded from `--config <toml/json>`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AppConfig {
    pub threads: Option<usize>,
    pub czti: CztiFileConfig,
    pub laxpc: LaxpcFileConfig,
    pub bary: BaryFileConfig,
    pub timing: TimingFileConfig,
    pub ui: UiFileConfig,
}

impl AppConfig {
    /// Validate enum fields and numeric bounds consistently across `astrosat --config`
    /// and `POST /api/config/parse`.
    pub fn validate(&self) -> Result<()> {
        if self.threads == Some(0) {
            return Err(CoreError::InvalidInput("threads must be >= 1".into()));
        }

        if let Some(es) = self.czti.entry_stage {
            if !(1..=3).contains(&es) {
                return Err(CoreError::InvalidInput(format!(
                    "czti.entry_stage must be 1..=3, got {es}"
                )));
            }
        }
        if let Some(xs) = self.czti.exit_stage {
            if !(1..=3).contains(&xs) {
                return Err(CoreError::InvalidInput(format!(
                    "czti.exit_stage must be 1..=3, got {xs}"
                )));
            }
        }
        if self.czti.entry_stage.unwrap_or(1) > self.czti.exit_stage.unwrap_or(3) {
            return Err(CoreError::InvalidInput(format!(
                "czti.entry_stage ({}) must be <= czti.exit_stage ({})",
                self.czti.entry_stage.unwrap_or(1),
                self.czti.exit_stage.unwrap_or(3)
            )));
        }
        if let Some(ref gt) = self.czti.gti_type {
            let norm = gt.trim().to_ascii_lowercase();
            if !matches!(norm.as_str(), "quad" | "common") {
                return Err(CoreError::InvalidInput(format!(
                    "invalid czti.gti_type '{}' (expected 'quad' or 'common')",
                    gt.trim()
                )));
            }
        }
        if let Some(ref al) = self.czti.algo {
            let norm = al.trim().to_ascii_lowercase();
            if !matches!(norm.as_str(), "default" | "generalized") {
                return Err(CoreError::InvalidInput(format!(
                    "invalid czti.algo '{}' (expected 'default' or 'generalized')",
                    al.trim()
                )));
            }
        }
        if let Some(ref ot) = self.czti.output_type {
            let norm = ot.trim().to_ascii_lowercase();
            if !matches!(
                norm.as_str(),
                "spec" | "lc" | "both" | "spectrum" | "lightcurve" | "all"
            ) {
                return Err(CoreError::InvalidInput(format!(
                    "invalid czti.output_type '{}' (expected 'spec', 'lc', or 'both')",
                    ot.trim()
                )));
            }
        }
        if let Some(ref ia) = self.czti.image_algo {
            let norm = ia.trim().to_ascii_lowercase();
            if !matches!(norm.as_str(), "balcc" | "fft") {
                return Err(CoreError::InvalidInput(format!(
                    "invalid czti.image_algo '{}' (expected 'balcc' or 'fft')",
                    ia.trim()
                )));
            }
        }
        if let Some(ov) = self.czti.oversampling {
            if ov == 0 || ov > 64 {
                return Err(CoreError::InvalidInput(format!(
                    "czti.oversampling must be in 1..=64, got {ov}"
                )));
            }
        }
        if let Some(ref qs) = self.czti.quads {
            if qs.is_empty() {
                return Err(CoreError::InvalidInput(
                    "czti.quads list cannot be empty".into(),
                ));
            }
            if qs.iter().any(|&q| q > 3) {
                return Err(CoreError::InvalidInput(
                    "czti.quads must contain values in 0..=3".into(),
                ));
            }
        }
        if let Some(ra) = self.czti.ra_deg {
            if !ra.is_finite() || !(-360.0..=360.0).contains(&ra) {
                return Err(CoreError::InvalidInput(
                    "czti.ra_deg must be finite in [-360, 360]".into(),
                ));
            }
        }
        if let Some(dec) = self.czti.dec_deg {
            if !dec.is_finite() || !(-90.0..=90.0).contains(&dec) {
                return Err(CoreError::InvalidInput(
                    "czti.dec_deg must be finite in [-90, 90]".into(),
                ));
            }
        }
        if let Some(tb) = self.czti.time_binsize {
            if !tb.is_finite() || tb <= 0.0 {
                return Err(CoreError::InvalidInput(
                    "czti.time_binsize must be finite and > 0".into(),
                ));
            }
        }
        if let Some(emin) = self.czti.emin_kev {
            if !emin.is_finite() || emin < 0.0 {
                return Err(CoreError::InvalidInput(
                    "czti.emin_kev must be finite and >= 0".into(),
                ));
            }
        }
        if let Some(emax) = self.czti.emax_kev {
            if !emax.is_finite() || emax <= 0.0 {
                return Err(CoreError::InvalidInput(
                    "czti.emax_kev must be finite and > 0".into(),
                ));
            }
        }
        if let (Some(emin), Some(emax)) = (self.czti.emin_kev, self.czti.emax_kev) {
            if emin >= emax {
                return Err(CoreError::InvalidInput(format!(
                    "czti.emin_kev ({emin}) must be < czti.emax_kev ({emax})"
                )));
            }
        }

        if let Some(ref us) = self.laxpc.units {
            if us.is_empty() {
                return Err(CoreError::InvalidInput(
                    "laxpc.units list cannot be empty".into(),
                ));
            }
            for u in us {
                if !is_valid_laxpc_unit_token(u) {
                    return Err(CoreError::InvalidInput(format!(
                        "invalid laxpc.units entry '{}' (expected '10', '20', '30', 'lxp1', 'lxp2', 'lxp3', or 'all')",
                        u.trim()
                    )));
                }
            }
        }
        if let Some(tb) = self.laxpc.tbin {
            if !tb.is_finite() || tb <= 0.0 {
                return Err(CoreError::InvalidInput(
                    "laxpc.tbin must be finite and > 0".into(),
                ));
            }
        }
        if let Some(dt) = self.laxpc.deadtime_sec {
            if !dt.is_finite() || dt < 0.0 {
                return Err(CoreError::InvalidInput(
                    "laxpc.deadtime_sec must be finite and >= 0".into(),
                ));
            }
        }
        if let Some(ch) = self.laxpc.spectrum_channels {
            if !matches!(ch, 256 | 512 | 1024) {
                return Err(CoreError::InvalidInput(format!(
                    "laxpc.spectrum_channels must be 256, 512, or 1024, got {ch}"
                )));
            }
        }
        if let Some(ref a) = self.laxpc.anodes {
            validate_laxpc_anodes_str(a)?;
        }
        if let Some(ref er) = self.laxpc.energy_range {
            validate_laxpc_energy_range_str(er)?;
        }
        if let Some(pmax) = self.laxpc.pha_max {
            if pmax > 1023 {
                return Err(CoreError::InvalidInput(format!(
                    "laxpc.pha_max must be <= 1023, got {pmax}"
                )));
            }
        }
        if let (Some(pmin), Some(pmax)) = (self.laxpc.pha_min, self.laxpc.pha_max) {
            if pmin > pmax {
                return Err(CoreError::InvalidInput(format!(
                    "laxpc.pha_min ({pmin}) must be <= laxpc.pha_max ({pmax})"
                )));
            }
        }
        if let Some(elv) = self.laxpc.min_elv_deg {
            if !elv.is_finite() {
                return Err(CoreError::InvalidInput(
                    "laxpc.min_elv_deg must be finite".into(),
                ));
            }
        }
        if let Some(lon) = self.laxpc.saa_entry_lon_deg {
            if !lon.is_finite() {
                return Err(CoreError::InvalidInput(
                    "laxpc.saa_entry_lon_deg must be finite".into(),
                ));
            }
        }
        if let Some(lon) = self.laxpc.saa_exit_lon_deg {
            if !lon.is_finite() {
                return Err(CoreError::InvalidInput(
                    "laxpc.saa_exit_lon_deg must be finite".into(),
                ));
            }
        }
        if let Some(gap) = self.laxpc.merge_gap_sec {
            if !gap.is_finite() || gap < 0.0 {
                return Err(CoreError::InvalidInput(
                    "laxpc.merge_gap_sec must be finite and >= 0".into(),
                ));
            }
        }
        if let Some(dur) = self.laxpc.min_duration_sec {
            if !dur.is_finite() || dur < 0.0 {
                return Err(CoreError::InvalidInput(
                    "laxpc.min_duration_sec must be finite and >= 0".into(),
                ));
            }
        }
        if let Some(ch) = self.laxpc.background_peak_ch {
            if !ch.is_finite() || ch <= 0.0 {
                return Err(CoreError::InvalidInput(
                    "laxpc.background_peak_ch must be finite and > 0".into(),
                ));
            }
        }

        if let Some(ref fr) = self.bary.frame {
            let norm = fr.trim().to_ascii_lowercase();
            if !matches!(norm.as_str(), "icrs" | "fk5") {
                return Err(CoreError::InvalidInput(format!(
                    "invalid bary.frame '{}' (expected 'icrs' or 'fk5')",
                    fr.trim()
                )));
            }
        }
        if let Some(ra) = self.bary.ra_deg {
            if !ra.is_finite() || !(-360.0..=360.0).contains(&ra) {
                return Err(CoreError::InvalidInput(
                    "bary.ra_deg must be finite in [-360, 360]".into(),
                ));
            }
        }
        if let Some(dec) = self.bary.dec_deg {
            if !dec.is_finite() || !(-90.0..=90.0).contains(&dec) {
                return Err(CoreError::InvalidInput(
                    "bary.dec_deg must be finite in [-90, 90]".into(),
                ));
            }
        }

        if let Some(h) = self.timing.harmonics {
            if h == 0 || h > 1024 {
                return Err(CoreError::InvalidInput(format!(
                    "timing.harmonics must be in 1..=1024, got {h}"
                )));
            }
        }
        if let Some(n) = self.timing.nbins {
            if !(2..=65_536).contains(&n) {
                return Err(CoreError::InvalidInput(format!(
                    "timing.nbins must be in 2..=65536, got {n}"
                )));
            }
        }

        if let Some(ref host) = self.ui.host {
            if host.trim().is_empty() {
                return Err(CoreError::InvalidInput(
                    "ui.host must be non-empty".into(),
                ));
            }
        }
        if let Some(port) = self.ui.port {
            if port == 0 {
                return Err(CoreError::InvalidInput(
                    "ui.port must be >= 1".into(),
                ));
            }
        }

        Ok(())
    }
}

fn is_valid_laxpc_unit_token(raw: &str) -> bool {
    matches!(
        raw.trim().to_ascii_lowercase().as_str(),
        "all"
            | "1"
            | "10"
            | "lxp1"
            | "lxp10"
            | "lx10"
            | "laxpc10"
            | "2"
            | "20"
            | "lxp2"
            | "lxp20"
            | "lx20"
            | "laxpc20"
            | "3"
            | "30"
            | "lxp3"
            | "lxp30"
            | "lx30"
            | "laxpc30"
    )
}

fn validate_laxpc_anodes_str(raw: &str) -> Result<()> {
    let t = raw.trim().to_ascii_lowercase();
    if matches!(
        t.as_str(),
        "all" | "a1-a7" | "layer1" | "l1" | "layer2" | "l2"
    ) {
        return Ok(());
    }
    if t.is_empty() {
        return Err(CoreError::InvalidInput(
            "laxpc.anodes cannot be empty".into(),
        ));
    }
    for part in t.split(',') {
        let p = part.trim().trim_start_matches('a');
        let id: u8 = p.parse().map_err(|_| {
            CoreError::InvalidInput(format!("invalid anode ID '{part}' in laxpc.anodes"))
        })?;
        if !(1..=10).contains(&id) {
            return Err(CoreError::InvalidInput(format!(
                "laxpc.anodes ID must be in 1..=10, got {id}"
            )));
        }
    }
    Ok(())
}

fn validate_laxpc_energy_range_str(raw: &str) -> Result<()> {
    let trimmed = raw.trim();
    let (elo, ehi) = if let Some((a, b)) = trimmed.split_once(':') {
        (a.trim(), b.trim())
    } else {
        let parts: Vec<&str> = trimmed.split('-').collect();
        if parts.len() != 2 {
            return Err(CoreError::InvalidInput(format!(
                "invalid laxpc.energy_range '{trimmed}' (expected EMIN:EMAX or EMIN-EMAX)"
            )));
        }
        (parts[0].trim(), parts[1].trim())
    };
    let lo: f64 = elo.parse().map_err(|_| {
        CoreError::InvalidInput(format!(
            "invalid lower bound '{elo}' in laxpc.energy_range"
        ))
    })?;
    let hi: f64 = ehi.parse().map_err(|_| {
        CoreError::InvalidInput(format!(
            "invalid upper bound '{ehi}' in laxpc.energy_range"
        ))
    })?;
    if !lo.is_finite() || !hi.is_finite() || lo < 0.0 || hi <= lo {
        return Err(CoreError::InvalidInput(
            "laxpc.energy_range must have finite 0 <= emin < emax".into(),
        ));
    }
    Ok(())
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CztiFileConfig {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_dir: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub l1_m0_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub l1_ss_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tct_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mkf_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub caldb_dir: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_dir: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_prefix: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mkf_threshold_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub user_gti_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub entry_stage: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub exit_stage: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gti_type: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub algo: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ra_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dec_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time_binsize: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub emin_kev: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub emax_kev: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_type: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mask_weight: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gen_resp: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub generate_eventfile: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub image_algo: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub oversampling: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quads: Option<Vec<usize>>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LaxpcFileConfig {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_dir: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_dir: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub caldb_dir: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub units: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tbin: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deadtime_sec: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub spectrum_channels: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub write_format_a_combined: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub write_format_b_events: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub write_orbit_file: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub write_ascii_products: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub single_only: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub layer1_only: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub anodes: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub energy_range: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pha_min: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pha_max: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_elv_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub saa_entry_lon_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub saa_exit_lon_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub merge_gap_sec: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_duration_sec: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub user_gti_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub user_gti_only: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub background_model_path: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub background_pha_path: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub background_peak_ch: Option<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct BaryFileConfig {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orbit_file: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ra_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dec_deg: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub frame: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ephemeris: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub write_checksums: Option<bool>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct TimingFileConfig {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub freq_hz: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fdot: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fddot: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub harmonics: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nbins: Option<usize>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct UiFileConfig {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub host: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub port: Option<u16>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub data_dir: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub open_browser: Option<bool>,
}

/// Quote and escape a string as a TOML basic string.
pub fn toml_quote(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Serialize an [`AppConfig`] to canonical TOML that round-trips identically via [`parse_toml_to_json`].
pub fn to_canonical_toml(cfg: &AppConfig) -> String {
    let mut out = String::new();
    if let Some(t) = cfg.threads {
        out.push_str(&format!("threads = {t}\n"));
    }

    fn push_header(out: &mut String, name: &str) {
        if !out.is_empty() && !out.ends_with("\n\n") {
            out.push('\n');
        }
        out.push_str(&format!("[{name}]\n"));
    }

    fn push_path(lines: &mut Vec<String>, key: &str, v: &Option<PathBuf>) {
        if let Some(p) = v {
            lines.push(format!("{key} = {}", toml_quote(&p.to_string_lossy())));
        }
    }
    fn push_str(lines: &mut Vec<String>, key: &str, v: &Option<String>) {
        if let Some(s) = v {
            lines.push(format!("{key} = {}", toml_quote(s)));
        }
    }
    fn push_disp<T: std::fmt::Display>(lines: &mut Vec<String>, key: &str, v: Option<T>) {
        if let Some(x) = v {
            lines.push(format!("{key} = {x}"));
        }
    }

    // [czti]
    {
        let c = &cfg.czti;
        let mut lines = Vec::new();
        push_path(&mut lines, "input_dir", &c.input_dir);
        push_path(&mut lines, "l1_m0_file", &c.l1_m0_file);
        push_path(&mut lines, "l1_ss_file", &c.l1_ss_file);
        push_path(&mut lines, "tct_file", &c.tct_file);
        push_path(&mut lines, "mkf_file", &c.mkf_file);
        push_path(&mut lines, "caldb_dir", &c.caldb_dir);
        push_path(&mut lines, "output_dir", &c.output_dir);
        push_str(&mut lines, "output_prefix", &c.output_prefix);
        push_path(&mut lines, "mkf_threshold_file", &c.mkf_threshold_file);
        push_path(&mut lines, "user_gti_file", &c.user_gti_file);
        push_disp(&mut lines, "entry_stage", c.entry_stage);
        push_disp(&mut lines, "exit_stage", c.exit_stage);
        push_str(&mut lines, "gti_type", &c.gti_type);
        push_str(&mut lines, "algo", &c.algo);
        push_disp(&mut lines, "ra_deg", c.ra_deg.filter(|x| x.is_finite()));
        push_disp(&mut lines, "dec_deg", c.dec_deg.filter(|x| x.is_finite()));
        push_disp(
            &mut lines,
            "time_binsize",
            c.time_binsize.filter(|x| x.is_finite()),
        );
        push_disp(&mut lines, "emin_kev", c.emin_kev.filter(|x| x.is_finite()));
        push_disp(&mut lines, "emax_kev", c.emax_kev.filter(|x| x.is_finite()));
        push_str(&mut lines, "output_type", &c.output_type);
        push_disp(&mut lines, "mask_weight", c.mask_weight);
        push_disp(&mut lines, "gen_resp", c.gen_resp);
        push_disp(&mut lines, "generate_eventfile", c.generate_eventfile);
        push_str(&mut lines, "image_algo", &c.image_algo);
        push_disp(&mut lines, "oversampling", c.oversampling);
        if let Some(ref qs) = c.quads {
            let arr = qs
                .iter()
                .map(|q| q.to_string())
                .collect::<Vec<_>>()
                .join(", ");
            lines.push(format!("quads = [{arr}]"));
        }
        if !lines.is_empty() {
            push_header(&mut out, "czti");
            for l in lines {
                out.push_str(&l);
                out.push('\n');
            }
        }
    }

    // [laxpc]
    {
        let l = &cfg.laxpc;
        let mut lines = Vec::new();
        push_path(&mut lines, "input_dir", &l.input_dir);
        push_path(&mut lines, "output_dir", &l.output_dir);
        push_path(&mut lines, "caldb_dir", &l.caldb_dir);
        if let Some(ref us) = l.units {
            let arr = us
                .iter()
                .map(|u| toml_quote(u))
                .collect::<Vec<_>>()
                .join(", ");
            lines.push(format!("units = [{arr}]"));
        }
        push_disp(&mut lines, "tbin", l.tbin.filter(|x| x.is_finite()));
        push_disp(
            &mut lines,
            "deadtime_sec",
            l.deadtime_sec.filter(|x| x.is_finite()),
        );
        push_disp(&mut lines, "spectrum_channels", l.spectrum_channels);
        push_disp(
            &mut lines,
            "write_format_a_combined",
            l.write_format_a_combined,
        );
        push_disp(&mut lines, "write_format_b_events", l.write_format_b_events);
        push_disp(&mut lines, "write_orbit_file", l.write_orbit_file);
        push_disp(&mut lines, "write_ascii_products", l.write_ascii_products);
        push_disp(&mut lines, "single_only", l.single_only);
        push_disp(&mut lines, "layer1_only", l.layer1_only);
        push_str(&mut lines, "anodes", &l.anodes);
        push_str(&mut lines, "energy_range", &l.energy_range);
        push_disp(&mut lines, "pha_min", l.pha_min);
        push_disp(&mut lines, "pha_max", l.pha_max);
        push_disp(
            &mut lines,
            "min_elv_deg",
            l.min_elv_deg.filter(|x| x.is_finite()),
        );
        push_disp(
            &mut lines,
            "saa_entry_lon_deg",
            l.saa_entry_lon_deg.filter(|x| x.is_finite()),
        );
        push_disp(
            &mut lines,
            "saa_exit_lon_deg",
            l.saa_exit_lon_deg.filter(|x| x.is_finite()),
        );
        push_disp(
            &mut lines,
            "merge_gap_sec",
            l.merge_gap_sec.filter(|x| x.is_finite()),
        );
        push_disp(
            &mut lines,
            "min_duration_sec",
            l.min_duration_sec.filter(|x| x.is_finite()),
        );
        push_path(&mut lines, "user_gti_file", &l.user_gti_file);
        push_disp(&mut lines, "user_gti_only", l.user_gti_only);
        push_path(&mut lines, "background_model_path", &l.background_model_path);
        push_path(&mut lines, "background_pha_path", &l.background_pha_path);
        push_disp(
            &mut lines,
            "background_peak_ch",
            l.background_peak_ch.filter(|x| x.is_finite()),
        );
        if !lines.is_empty() {
            push_header(&mut out, "laxpc");
            for line in lines {
                out.push_str(&line);
                out.push('\n');
            }
        }
    }

    // [bary]
    {
        let b = &cfg.bary;
        let mut lines = Vec::new();
        push_path(&mut lines, "orbit_file", &b.orbit_file);
        push_disp(&mut lines, "ra_deg", b.ra_deg.filter(|x| x.is_finite()));
        push_disp(&mut lines, "dec_deg", b.dec_deg.filter(|x| x.is_finite()));
        push_str(&mut lines, "frame", &b.frame);
        push_str(&mut lines, "ephemeris", &b.ephemeris);
        push_disp(&mut lines, "write_checksums", b.write_checksums);
        if !lines.is_empty() {
            push_header(&mut out, "bary");
            for line in lines {
                out.push_str(&line);
                out.push('\n');
            }
        }
    }

    // [timing]
    {
        let t = &cfg.timing;
        let mut lines = Vec::new();
        push_disp(&mut lines, "freq_hz", t.freq_hz.filter(|x| x.is_finite()));
        push_disp(&mut lines, "fdot", t.fdot.filter(|x| x.is_finite()));
        push_disp(&mut lines, "fddot", t.fddot.filter(|x| x.is_finite()));
        push_disp(&mut lines, "harmonics", t.harmonics);
        push_disp(&mut lines, "nbins", t.nbins);
        if !lines.is_empty() {
            push_header(&mut out, "timing");
            for line in lines {
                out.push_str(&line);
                out.push('\n');
            }
        }
    }

    // [ui]
    {
        let u = &cfg.ui;
        let mut lines = Vec::new();
        push_str(&mut lines, "host", &u.host);
        push_disp(&mut lines, "port", u.port);
        push_path(&mut lines, "data_dir", &u.data_dir);
        push_disp(&mut lines, "open_browser", u.open_browser);
        if !lines.is_empty() {
            push_header(&mut out, "ui");
            for line in lines {
                out.push_str(&line);
                out.push('\n');
            }
        }
    }

    out
}

/// Parse a TOML document into a `serde_json::Value::Object`.
pub fn parse_toml_to_json(input: &str) -> Result<Value> {
    let input = input.strip_prefix('\u{FEFF}').unwrap_or(input);
    let mut root = Map::new();
    let mut current_table: Vec<String> = Vec::new();
    let mut pending_array: Option<(usize, Vec<String>, String)> = None;

    for (line_idx, raw_line) in input.lines().enumerate() {
        let line = strip_toml_comment(raw_line).trim();

        if let Some((start_idx, full_key, mut buf)) = pending_array.take() {
            if !line.is_empty() {
                if !buf.is_empty() {
                    buf.push(' ');
                }
                buf.push_str(line);
            }
            match array_bracket_balance(&buf) {
                Ok(0) => {
                    let trimmed = buf.trim();
                    let val = parse_toml_value(trimmed).map_err(|e| {
                        CoreError::InvalidInput(format!(
                            "line {}: invalid value '{}': {e}",
                            start_idx + 1,
                            trimmed
                        ))
                    })?;
                    insert_nested(&mut root, &full_key, val)?;
                }
                Ok(_) => {
                    pending_array = Some((start_idx, full_key, buf));
                }
                Err(e) => {
                    return Err(CoreError::InvalidInput(format!(
                        "line {}: {e}",
                        line_idx + 1
                    )));
                }
            }
            continue;
        }

        if line.is_empty() {
            continue;
        }

        if line.starts_with('[') && line.ends_with(']') {
            let section = line[1..line.len() - 1].trim();
            if section.is_empty() || section.starts_with('[') {
                return Err(CoreError::InvalidInput(format!(
                    "line {}: unsupported TOML header '{}'",
                    line_idx + 1,
                    line
                )));
            }
            let mut segs = Vec::new();
            for s in section.split('.') {
                let seg = s.trim().trim_matches('"').trim_matches('\'');
                if seg.is_empty() {
                    return Err(CoreError::InvalidInput(format!(
                        "line {}: empty section segment in '{}'",
                        line_idx + 1,
                        line
                    )));
                }
                segs.push(seg.to_string());
            }
            current_table = segs;
            ensure_table_path(&mut root, &current_table)?;
            continue;
        }

        let Some((lhs, rhs)) = split_key_value(line) else {
            return Err(CoreError::InvalidInput(format!(
                "line {}: expected key = value, got '{}'",
                line_idx + 1,
                line
            )));
        };
        let mut full_key: Vec<String> = current_table.clone();
        for part in lhs.split('.') {
            let k = part.trim().trim_matches('"').trim_matches('\'');
            if k.is_empty() {
                return Err(CoreError::InvalidInput(format!(
                    "line {}: empty key segment",
                    line_idx + 1
                )));
            }
            full_key.push(k.to_string());
        }
        let rhs_trimmed = rhs.trim();
        if rhs_trimmed.starts_with('[') {
            match array_bracket_balance(rhs_trimmed) {
                Ok(d) if d > 0 => {
                    pending_array = Some((line_idx, full_key, rhs_trimmed.to_string()));
                    continue;
                }
                Ok(_) => {}
                Err(e) => {
                    return Err(CoreError::InvalidInput(format!(
                        "line {}: invalid value '{}': {e}",
                        line_idx + 1,
                        rhs_trimmed
                    )));
                }
            }
        }
        let val = parse_toml_value(rhs_trimmed).map_err(|e| {
            CoreError::InvalidInput(format!(
                "line {}: invalid value '{}': {e}",
                line_idx + 1,
                rhs_trimmed
            ))
        })?;
        insert_nested(&mut root, &full_key, val)?;
    }

    if let Some((start_idx, _, _)) = pending_array {
        return Err(CoreError::InvalidInput(format!(
            "line {}: unterminated TOML array",
            start_idx + 1
        )));
    }

    Ok(Value::Object(root))
}

fn array_bracket_balance(s: &str) -> Result<usize> {
    let mut depth = 0usize;
    let mut in_dq = false;
    let mut in_sq = false;
    let mut escape = false;
    for ch in s.chars() {
        if escape {
            escape = false;
            continue;
        }
        if in_dq {
            if ch == '\\' {
                escape = true;
            } else if ch == '"' {
                in_dq = false;
            }
        } else if in_sq {
            if ch == '\'' {
                in_sq = false;
            }
        } else {
            match ch {
                '"' => in_dq = true,
                '\'' => in_sq = true,
                '[' => depth += 1,
                ']' => {
                    if depth == 0 {
                        return Err(CoreError::InvalidInput(
                            "unmatched ']' in TOML array".into(),
                        ));
                    }
                    depth -= 1;
                }
                _ => {}
            }
        }
    }
    if in_dq || in_sq {
        return Err(CoreError::InvalidInput(
            "unterminated string in TOML array".into(),
        ));
    }
    Ok(depth)
}

fn strip_toml_comment(line: &str) -> &str {
    let mut in_dq = false;
    let mut in_sq = false;
    let mut escape = false;
    for (i, ch) in line.char_indices() {
        if escape {
            escape = false;
            continue;
        }
        if in_dq {
            if ch == '\\' {
                escape = true;
            } else if ch == '"' {
                in_dq = false;
            }
        } else if in_sq {
            if ch == '\'' {
                in_sq = false;
            }
        } else {
            match ch {
                '"' => in_dq = true,
                '\'' => in_sq = true,
                '#' => return &line[..i],
                _ => {}
            }
        }
    }
    line
}

fn split_key_value(line: &str) -> Option<(&str, &str)> {
    let mut in_dq = false;
    let mut in_sq = false;
    let mut escape = false;
    for (i, ch) in line.char_indices() {
        if escape {
            escape = false;
            continue;
        }
        if in_dq {
            if ch == '\\' {
                escape = true;
            } else if ch == '"' {
                in_dq = false;
            }
        } else if in_sq {
            if ch == '\'' {
                in_sq = false;
            }
        } else if ch == '=' {
            return Some((&line[..i], &line[i + 1..]));
        } else if ch == '"' {
            in_dq = true;
        } else if ch == '\'' {
            in_sq = true;
        }
    }
    None
}

fn parse_toml_basic_string(s: &str) -> Result<String> {
    let mut chars = s.chars();
    if chars.next() != Some('"') {
        return Err(CoreError::InvalidInput("expected opening quote".into()));
    }
    let mut out = String::with_capacity(s.len().saturating_sub(2));
    let mut closed = false;
    while let Some(ch) = chars.next() {
        match ch {
            '"' => {
                closed = true;
                break;
            }
            '\\' => match chars.next() {
                Some('"') => out.push('"'),
                Some('\\') => out.push('\\'),
                Some('n') => out.push('\n'),
                Some('t') => out.push('\t'),
                Some('r') => out.push('\r'),
                Some(other) => {
                    return Err(CoreError::InvalidInput(format!(
                        "unsupported escape sequence '\\{other}' in TOML string"
                    )));
                }
                None => {
                    return Err(CoreError::InvalidInput(
                        "unterminated escape sequence in TOML string".into(),
                    ));
                }
            },
            _ => out.push(ch),
        }
    }
    if !closed || chars.next().is_some() {
        return Err(CoreError::InvalidInput(format!(
            "malformed TOML string '{s}'"
        )));
    }
    Ok(out)
}

fn parse_toml_value(s: &str) -> Result<Value> {
    if s.starts_with('"') {
        return Ok(Value::String(parse_toml_basic_string(s)?));
    }
    if s.starts_with('\'')
        && s.ends_with('\'')
        && s.len() >= 2
        && !s[1..s.len() - 1].contains('\'')
    {
        return Ok(Value::String(s[1..s.len() - 1].to_string()));
    }
    if s == "true" {
        return Ok(Value::Bool(true));
    }
    if s == "false" {
        return Ok(Value::Bool(false));
    }
    if s.starts_with('[') && s.ends_with(']') {
        let inner = s[1..s.len() - 1].trim();
        if inner.is_empty() {
            return Ok(Value::Array(Vec::new()));
        }
        let items = split_toml_array(inner)?;
        let mut arr = Vec::with_capacity(items.len());
        for (idx, item) in items.iter().enumerate() {
            let trimmed = item.trim();
            if trimmed.is_empty() {
                if idx + 1 == items.len() && !arr.is_empty() {
                    continue;
                }
                return Err(CoreError::InvalidInput(
                    "empty element in TOML array".into(),
                ));
            }
            arr.push(parse_toml_value(trimmed)?);
        }
        return Ok(Value::Array(arr));
    }
    let clean_num = s.replace('_', "");
    if let Ok(i) = clean_num.parse::<i64>() {
        return Ok(Value::Number(i.into()));
    }
    if let Ok(f) = clean_num.parse::<f64>() {
        let num = serde_json::Number::from_f64(f)
            .ok_or_else(|| CoreError::InvalidInput(format!("non-finite float '{s}'")))?;
        return Ok(Value::Number(num));
    }
    Err(CoreError::InvalidInput(format!(
        "unrecognized TOML value '{s}'"
    )))
}

fn split_toml_array(inner: &str) -> Result<Vec<&str>> {
    let mut parts = Vec::new();
    let mut start = 0usize;
    let mut depth = 0usize;
    let mut in_dq = false;
    let mut in_sq = false;
    let mut escape = false;

    for (i, ch) in inner.char_indices() {
        if escape {
            escape = false;
            continue;
        }
        if in_dq {
            if ch == '\\' {
                escape = true;
            } else if ch == '"' {
                in_dq = false;
            }
        } else if in_sq {
            if ch == '\'' {
                in_sq = false;
            }
        } else {
            match ch {
                '"' => in_dq = true,
                '\'' => in_sq = true,
                '[' => depth += 1,
                ']' => {
                    if depth == 0 {
                        return Err(CoreError::InvalidInput(
                            "unmatched ']' in TOML array".into(),
                        ));
                    }
                    depth -= 1;
                }
                ',' if depth == 0 => {
                    parts.push(&inner[start..i]);
                    start = i + 1;
                }
                _ => {}
            }
        }
    }
    if in_dq || in_sq || depth != 0 {
        return Err(CoreError::InvalidInput(
            "unterminated string or bracket in TOML array".into(),
        ));
    }
    parts.push(&inner[start..]);
    Ok(parts)
}

fn ensure_table_path(root: &mut Map<String, Value>, path: &[String]) -> Result<()> {
    let mut cur = root;
    for seg in path {
        let entry = cur
            .entry(seg.clone())
            .or_insert_with(|| Value::Object(Map::new()));
        cur = entry.as_object_mut().ok_or_else(|| {
            CoreError::InvalidInput(format!("key '{seg}' conflicts with non-table value"))
        })?;
    }
    Ok(())
}

fn insert_nested(root: &mut Map<String, Value>, path: &[String], val: Value) -> Result<()> {
    if path.is_empty() {
        return Err(CoreError::InvalidInput("empty key path".into()));
    }
    let mut cur = root;
    for seg in &path[..path.len() - 1] {
        let entry = cur
            .entry(seg.clone())
            .or_insert_with(|| Value::Object(Map::new()));
        cur = entry.as_object_mut().ok_or_else(|| {
            CoreError::InvalidInput(format!("key '{seg}' conflicts with non-table value"))
        })?;
    }
    cur.insert(path[path.len() - 1].clone(), val);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_sections_dotted_keys_arrays_escapes_and_comments() {
        let doc = r#"
# Comment
threads = 4
timing.fddot = 1.5e-20
timing.nbins = 1_024

[czti]
output_prefix = "crab\n\"q0\"\t\\ok\r" # trailing comment
quads = [0, 1, 2, 3]
oversampling = 32
mask_weight = true

[laxpc]
units = ['10', '20']
energy_range = "3:80" # inline # hash inside quotes: "a#b"
background_model_path = "C:\\data\\backfit.lxp10"
"#;
        let val = parse_toml_to_json(doc).unwrap();
        let cfg: AppConfig = serde_json::from_value(val).unwrap();
        assert_eq!(cfg.threads, Some(4));
        assert_eq!(cfg.timing.nbins, Some(1024));
        assert!((cfg.timing.fddot.unwrap() - 1.5e-20).abs() < 1e-30);
        assert_eq!(cfg.czti.output_prefix.as_deref(), Some("crab\n\"q0\"\t\\ok\r"));
        assert_eq!(cfg.czti.quads.as_deref(), Some(&[0, 1, 2, 3][..]));
        assert_eq!(cfg.czti.oversampling, Some(32));
        assert_eq!(cfg.czti.mask_weight, Some(true));
        assert_eq!(
            cfg.laxpc.units.as_deref(),
            Some(&["10".to_string(), "20".to_string()][..])
        );
        assert_eq!(cfg.laxpc.energy_range.as_deref(), Some("3:80"));
        assert_eq!(
            cfg.laxpc.background_model_path.as_deref(),
            Some(std::path::Path::new(r"C:\data\backfit.lxp10"))
        );

        let canon = to_canonical_toml(&cfg);
        let reparsed: AppConfig =
            serde_json::from_value(parse_toml_to_json(&canon).unwrap()).unwrap();
        assert_eq!(reparsed.threads, cfg.threads);
        assert_eq!(reparsed.czti.output_prefix, cfg.czti.output_prefix);
        assert_eq!(reparsed.czti.oversampling, cfg.czti.oversampling);
        assert_eq!(reparsed.laxpc.units, cfg.laxpc.units);
        assert_eq!(reparsed.timing.nbins, cfg.timing.nbins);
    }

    #[test]
    fn rejects_malformed_toml_and_unknown_schema_fields() {
        let bad_cases = [
            "[]",
            "[[array_table]]",
            "[czti.]",
            "not_a_kv",
            "[czti]\n= 1",
            "[czti]\ns = \"bad\\x\"",
            "[czti]\ns = \"unterminated",
            "[czti]\na = [1, 2",
            "[czti]\na = [\"unterminated]",
            "[czti]\nf = nan",
            "a = 1\na.b = 2",
        ];
        for doc in bad_cases {
            assert!(
                parse_toml_to_json(doc).is_err(),
                "expected TOML error for: {doc}"
            );
        }

        // Unknown field and unknown section rejected by AppConfig
        let unknown_key = parse_toml_to_json("[czti]\nbogus_key = 42").unwrap();
        assert!(serde_json::from_value::<AppConfig>(unknown_key).is_err());

        let unknown_sec = parse_toml_to_json("[bogus]\nx = 1").unwrap();
        assert!(serde_json::from_value::<AppConfig>(unknown_sec).is_err());

        // Consecutive empty commas inside arrays must fail
        for bad_arr in [
            "[czti]\nquads = [0, , 1]",
            "[czti]\nquads = [, 0]",
            "[czti]\nquads = [,]",
            "[czti]\nquads = [0, ,]",
        ] {
            assert!(
                parse_toml_to_json(bad_arr).is_err(),
                "expected array error for: {bad_arr}"
            );
        }
    }

    #[test]
    fn supports_utf8_bom_multiline_arrays_with_comments_and_enum_validation() {
        let doc = "\u{FEFF}threads = 2\n[czti]\nquads = [\n  # first\n  0, # q0\n  1,\n  2, # trailing comma\n]\n[bary]\nframe = \"fk5\"\n";
        let val = parse_toml_to_json(doc).unwrap();
        let cfg: AppConfig = serde_json::from_value(val).unwrap();
        cfg.validate().unwrap();
        assert_eq!(cfg.threads, Some(2));
        assert_eq!(cfg.czti.quads.as_deref(), Some(&[0, 1, 2][..]));
        assert_eq!(cfg.bary.frame.as_deref(), Some("fk5"));

        // Trailing comma on single-line array still works
        let trailing = parse_toml_to_json("[czti]\nquads = [0, 1,]").unwrap();
        let trailing_cfg: AppConfig = serde_json::from_value(trailing).unwrap();
        assert_eq!(trailing_cfg.czti.quads.as_deref(), Some(&[0, 1][..]));

        // Invalid enums rejected consistently by AppConfig::validate()
        for bad_cfg_toml in [
            "[bary]\nframe = \"galactic\"",
            "[czti]\nalgo = \"bogus\"",
            "[czti]\ngti_type = \"invalid\"",
            "[czti]\noutput_type = \"bogus\"",
            "[czti]\nimage_algo = \"psf\"",
            "[czti]\nentry_stage = 3\nexit_stage = 1",
            "[laxpc]\nspectrum_channels = 128",
        ] {
            let v = parse_toml_to_json(bad_cfg_toml).unwrap();
            let c: AppConfig = serde_json::from_value(v).unwrap();
            assert!(
                c.validate().is_err(),
                "expected AppConfig::validate() error for:\n{bad_cfg_toml}"
            );
        }
    }

    #[test]
    fn validates_laxpc_czti_and_ui_bounds_consistently() {
        let valid_toml = r#"
[czti]
ra_deg = 83.6331
dec_deg = -29.0078
time_binsize = 0.5

[laxpc]
units = ["10", "lxp2", "all"]
tbin = 0.25
deadtime_sec = 0.0
anodes = "a1, a2, layer1"
energy_range = "3-20"
merge_gap_sec = 0.0
min_duration_sec = 0.0
background_peak_ch = 213.5

[ui]
host = "127.0.0.1"
port = 8085
"#;
        // Note: "a1, a2, layer1" isn't a valid anode ID list; test valid forms separately
        let _ = valid_toml;
        for good_cfg in [
            "[laxpc]\nenergy_range = \"3:20\"\nunits = [\"10\", \"lxp2\", \"30\"]\ntbin = 0.5\ndeadtime_sec = 0.0\nanodes = \"all\"\nmerge_gap_sec = 0.0\nmin_duration_sec = 0.0\nbackground_peak_ch = 213.0\n",
            "[laxpc]\nenergy_range = \"0-80.5\"\nunits = [\"all\"]\nanodes = \"1,2,3\"\n",
            "[laxpc]\nanodes = \"layer1\"\n",
            "[laxpc]\nanodes = \"a1, a7\"\n",
            "[czti]\nra_deg = 360.0\ndec_deg = -90.0\ntime_binsize = 1.0\n",
            "[ui]\nhost = \"0.0.0.0\"\nport = 8085\n",
        ] {
            let v = parse_toml_to_json(good_cfg).unwrap();
            let c: AppConfig = serde_json::from_value(v).unwrap();
            assert!(
                c.validate().is_ok(),
                "expected AppConfig::validate() OK for:\n{good_cfg}"
            );
        }

        for bad_cfg in [
            "[laxpc]\nenergy_range = \"bogus\"\n",
            "[laxpc]\nenergy_range = \"80:10\"\n",
            "[laxpc]\nenergy_range = \"20-3\"\n",
            "[laxpc]\nenergy_range = \"-5:20\"\n",
            "[laxpc]\nenergy_range = \"10:10\"\n",
            "[laxpc]\nenergy_range = \"\"\n",
            "[laxpc]\nunits = [\"99\"]\n",
            "[laxpc]\nunits = [\"\"]\n",
            "[laxpc]\nunits = []\n",
            "[laxpc]\ntbin = -1.0\n",
            "[laxpc]\ntbin = 0.0\n",
            "[laxpc]\ndeadtime_sec = -0.1\n",
            "[laxpc]\nanodes = \"bogus\"\n",
            "[laxpc]\nanodes = \"1,99\"\n",
            "[laxpc]\nanodes = \"0\"\n",
            "[laxpc]\nanodes = \"\"\n",
            "[laxpc]\nanodes = \"1,,2\"\n",
            "[laxpc]\nmerge_gap_sec = -1.0\n",
            "[laxpc]\nmin_duration_sec = -0.5\n",
            "[laxpc]\nbackground_peak_ch = 0.0\n",
            "[laxpc]\nbackground_peak_ch = -10.0\n",
            "[czti]\ntime_binsize = -1.0\n",
            "[czti]\ntime_binsize = 0.0\n",
            "[czti]\nra_deg = 400.0\n",
            "[czti]\nra_deg = -400.0\n",
            "[czti]\ndec_deg = 120.0\n",
            "[czti]\ndec_deg = -91.0\n",
            "[czti]\nquads = []\n",
            "[ui]\nport = 0\n",
            "[ui]\nhost = \"\"\n",
            "[ui]\nhost = \"   \"\n",
        ] {
            let v = parse_toml_to_json(bad_cfg).unwrap();
            let c: AppConfig = serde_json::from_value(v).unwrap();
            assert!(
                c.validate().is_err(),
                "expected AppConfig::validate() error for:\n{bad_cfg}"
            );
        }
    }

    #[test]
    fn rejects_empty_czti_quads_while_allowing_omitted() {
        let omitted_val = parse_toml_to_json("[czti]\nra_deg = 83.6331\n").unwrap();
        let omitted_cfg: AppConfig = serde_json::from_value(omitted_val).unwrap();
        assert!(omitted_cfg.czti.quads.is_none());
        omitted_cfg.validate().unwrap();

        let empty_val = parse_toml_to_json("[czti]\nquads = []\n").unwrap();
        let empty_cfg: AppConfig = serde_json::from_value(empty_val).unwrap();
        let err = empty_cfg.validate().unwrap_err();
        assert!(
            err.to_string().contains("czti.quads list cannot be empty"),
            "unexpected error for empty czti.quads: {err}"
        );
    }
}
