//! Configuration file loader supporting both JSON (`.json`) and TOML (`.toml`) formats.
//!
//! Uses a lightweight pure-Rust TOML subset parser (tables `[section]`, nested dotted keys,
//! comments `#`, strings, integers, floats, booleans, and inline arrays) so zero extra external
//! dependencies are required beyond `serde_json`.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

/// Top-level CLI configuration loaded from `--config <toml/json>`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct AppConfig {
    pub threads: Option<usize>,
    pub czti: CztiFileConfig,
    pub laxpc: LaxpcFileConfig,
    pub bary: BaryFileConfig,
    pub timing: TimingFileConfig,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct CztiFileConfig {
    pub input_dir: Option<PathBuf>,
    pub l1_m0_file: Option<PathBuf>,
    pub l1_ss_file: Option<PathBuf>,
    pub tct_file: Option<PathBuf>,
    pub mkf_file: Option<PathBuf>,
    pub caldb_dir: Option<PathBuf>,
    pub output_dir: Option<PathBuf>,
    pub output_prefix: Option<String>,
    pub mkf_threshold_file: Option<PathBuf>,
    pub user_gti_file: Option<PathBuf>,
    pub entry_stage: Option<u8>,
    pub exit_stage: Option<u8>,
    pub gti_type: Option<String>,
    pub algo: Option<String>,
    pub ra_deg: Option<f64>,
    pub dec_deg: Option<f64>,
    pub time_binsize: Option<f64>,
    pub emin_kev: Option<f32>,
    pub emax_kev: Option<f32>,
    pub output_type: Option<String>,
    pub mask_weight: Option<bool>,
    pub gen_resp: Option<bool>,
    pub generate_eventfile: Option<bool>,
    pub image_algo: Option<String>,
    pub oversampling: Option<usize>,
    pub quads: Option<Vec<usize>>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct LaxpcFileConfig {
    pub input_dir: Option<PathBuf>,
    pub output_dir: Option<PathBuf>,
    pub caldb_dir: Option<PathBuf>,
    pub units: Option<Vec<String>>,
    pub tbin: Option<f64>,
    pub deadtime_sec: Option<f64>,
    pub spectrum_channels: Option<usize>,
    pub write_format_a_combined: Option<bool>,
    pub write_format_b_events: Option<bool>,
    pub write_orbit_file: Option<bool>,
    pub write_ascii_products: Option<bool>,
    pub single_only: Option<bool>,
    pub layer1_only: Option<bool>,
    pub anodes: Option<String>,
    pub energy_range: Option<String>,
    pub pha_min: Option<u16>,
    pub pha_max: Option<u16>,
    pub min_elv_deg: Option<f64>,
    pub saa_entry_lon_deg: Option<f64>,
    pub saa_exit_lon_deg: Option<f64>,
    pub merge_gap_sec: Option<f64>,
    pub min_duration_sec: Option<f64>,
    pub user_gti_file: Option<PathBuf>,
    pub user_gti_only: Option<bool>,
    pub background_model_path: Option<PathBuf>,
    pub background_pha_path: Option<PathBuf>,
    pub background_peak_ch: Option<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct BaryFileConfig {
    pub orbit_file: Option<PathBuf>,
    pub ra_deg: Option<f64>,
    pub dec_deg: Option<f64>,
    pub frame: Option<String>,
    pub ephemeris: Option<String>,
    pub write_checksums: Option<bool>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct TimingFileConfig {
    pub freq_hz: Option<f64>,
    pub fdot: Option<f64>,
    pub fddot: Option<f64>,
    pub harmonics: Option<usize>,
    pub nbins: Option<usize>,
}

/// Load [`AppConfig`] from a `.json` or `.toml` file.
pub fn load_app_config(path: &Path) -> Result<AppConfig> {
    let raw = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read config file {}", path.display()))?;
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    let val = if ext == "json" || raw.trim_start().starts_with('{') {
        serde_json::from_str::<Value>(&raw)
            .with_context(|| format!("invalid JSON in config file {}", path.display()))?
    } else {
        parse_toml_to_json(&raw)
            .with_context(|| format!("invalid TOML in config file {}", path.display()))?
    };

    serde_json::from_value::<AppConfig>(val)
        .with_context(|| format!("config schema mismatch in {}", path.display()))
}

/// Parse a TOML document into a `serde_json::Value::Object`.
pub fn parse_toml_to_json(input: &str) -> Result<Value> {
    let mut root = Map::new();
    let mut current_table: Vec<String> = Vec::new();

    for (line_idx, raw_line) in input.lines().enumerate() {
        let line = strip_toml_comment(raw_line).trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('[') && line.ends_with(']') {
            let section = line[1..line.len() - 1].trim();
            if section.is_empty() || section.starts_with('[') {
                bail!("line {}: unsupported TOML header '{}'", line_idx + 1, line);
            }
            current_table = section
                .split('.')
                .map(|s| s.trim().trim_matches('"').trim_matches('\'').to_string())
                .collect();
            ensure_table_path(&mut root, &current_table)?;
            continue;
        }

        let Some((lhs, rhs)) = split_key_value(line) else {
            bail!("line {}: expected key = value, got '{}'", line_idx + 1, line);
        };
        let mut full_key: Vec<String> = current_table.clone();
        for part in lhs.split('.') {
            let k = part.trim().trim_matches('"').trim_matches('\'');
            if k.is_empty() {
                bail!("line {}: empty key segment", line_idx + 1);
            }
            full_key.push(k.to_string());
        }
        let val = parse_toml_value(rhs.trim())
            .with_context(|| format!("line {}: invalid value '{}'", line_idx + 1, rhs.trim()))?;
        insert_nested(&mut root, &full_key, val)?;
    }

    Ok(Value::Object(root))
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
        bail!("expected opening quote");
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
                Some(other) => bail!("unsupported escape sequence '\\{other}' in TOML string"),
                None => bail!("unterminated escape sequence in TOML string"),
            },
            _ => out.push(ch),
        }
    }
    if !closed || chars.next().is_some() {
        bail!("malformed TOML string '{s}'");
    }
    Ok(out)
}

fn parse_toml_value(s: &str) -> Result<Value> {
    if s.starts_with('"') {
        return Ok(Value::String(parse_toml_basic_string(s)?));
    }
    if s.starts_with('\'') && s.ends_with('\'') && s.len() >= 2 && !s[1..s.len() - 1].contains('\'') {
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
        for item in items {
            let trimmed = item.trim();
            if !trimmed.is_empty() {
                arr.push(parse_toml_value(trimmed)?);
            }
        }
        return Ok(Value::Array(arr));
    }
    let clean_num = s.replace('_', "");
    if let Ok(i) = clean_num.parse::<i64>() {
        return Ok(Value::Number(i.into()));
    }
    if let Ok(f) = clean_num.parse::<f64>() {
        let num = serde_json::Number::from_f64(f)
            .with_context(|| format!("non-finite float '{s}'"))?;
        return Ok(Value::Number(num));
    }
    bail!("unrecognized TOML value '{s}'")
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
                ']' => depth = depth.saturating_sub(1),
                ',' if depth == 0 => {
                    parts.push(&inner[start..i]);
                    start = i + 1;
                }
                _ => {}
            }
        }
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
        cur = entry
            .as_object_mut()
            .with_context(|| format!("key '{seg}' conflicts with non-table value"))?;
    }
    Ok(())
}

fn insert_nested(root: &mut Map<String, Value>, path: &[String], val: Value) -> Result<()> {
    if path.is_empty() {
        bail!("empty key path");
    }
    let mut cur = root;
    for seg in &path[..path.len() - 1] {
        let entry = cur
            .entry(seg.clone())
            .or_insert_with(|| Value::Object(Map::new()));
        cur = entry
            .as_object_mut()
            .with_context(|| format!("key '{seg}' conflicts with non-table value"))?;
    }
    cur.insert(path[path.len() - 1].clone(), val);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_toml_sections_arrays_and_comments() {
        let doc = r#"
# Global
threads = 8

[laxpc]
input_dir = "/path/to/data/sample/laxpc_l1"
units = ["10", "20", "30"]
tbin = 0.5
write_format_a_combined = true

[czti]
quads = [0, 1, 2, 3]
ra_deg = 83.6331 # Crab RA
"#;
        let json = parse_toml_to_json(doc).unwrap();
        let cfg: AppConfig = serde_json::from_value(json).unwrap();
        assert_eq!(cfg.threads, Some(8));
        assert_eq!(
            cfg.laxpc.units.unwrap(),
            vec!["10".to_string(), "20".to_string(), "30".to_string()]
        );
        assert_eq!(cfg.laxpc.tbin, Some(0.5));
        assert_eq!(cfg.laxpc.write_format_a_combined, Some(true));
        assert_eq!(cfg.czti.quads.unwrap(), vec![0, 1, 2, 3]);
        assert!((cfg.czti.ra_deg.unwrap() - 83.6331).abs() < 1e-6);

        let escaped = r#"
[laxpc]
background_pha_path = "C:\\data\\back\"quoted\".pha"
background_peak_ch = 213.0
energy_range = "3:20"
"#;
        let esc_cfg: AppConfig =
            serde_json::from_value(parse_toml_to_json(escaped).unwrap()).unwrap();
        assert_eq!(
            esc_cfg.laxpc.background_pha_path.unwrap().to_str().unwrap(),
            r#"C:\data\back"quoted".pha"#
        );
        assert_eq!(esc_cfg.laxpc.background_peak_ch, Some(213.0));
        assert_eq!(esc_cfg.laxpc.energy_range.as_deref(), Some("3:20"));
    }
}
