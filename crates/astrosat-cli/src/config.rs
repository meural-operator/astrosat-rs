//! Configuration file loader supporting both JSON (`.json`) and TOML (`.toml`) formats.
//!
//! Delegates TOML parsing and [`AppConfig`] schema definitions to `astrosat_core::toml_lite`
//! so `astrosat --config` and `astrosat-ui` (`POST /api/config/parse`) share identical semantics.

use std::path::Path;

use anyhow::{Context, Result};
pub use astrosat_core::{
    parse_toml_to_json, AppConfig, BaryFileConfig, CztiFileConfig, LaxpcFileConfig,
};
use serde_json::Value;

/// Load [`AppConfig`] from a `.json` or `.toml` file.
pub fn load_app_config(path: &Path) -> Result<AppConfig> {
    let raw = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read config file {}", path.display()))?;
    let raw_trimmed = raw.strip_prefix('\u{FEFF}').unwrap_or(&raw);
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    let val = if ext == "json" || raw_trimmed.trim_start().starts_with('{') {
        serde_json::from_str::<Value>(raw_trimmed)
            .with_context(|| format!("invalid JSON in config file {}", path.display()))?
    } else {
        parse_toml_to_json(raw_trimmed)
            .with_context(|| format!("invalid TOML in config file {}", path.display()))?
    };

    let cfg = serde_json::from_value::<AppConfig>(val)
        .with_context(|| format!("config schema mismatch in {}", path.display()))?;
    cfg.validate()
        .with_context(|| format!("invalid config value in {}", path.display()))?;
    Ok(cfg)
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
input_dir = "/mnt/d/ISRO/data/sample/laxpc_l1"
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
