//! Zero-external-runtime-dependency local web UI and scientific REST API for `astrosat-rs`.
//!
//! Serves an embedded offline HTML5/CSS3/JS scientific workbench alongside multi-core `rayon`
//! endpoints for FITS inspection & OGIP validation, light curves, spectra, 2D coded-aperture
//! imaging + WCS, multi-band phase folding & pulse alignment, and in-process pipeline execution.

pub mod browse;
pub mod config_api;
pub mod fits_api;
pub mod http_util;
pub mod jobs;
pub mod pulse;
pub mod science;
pub mod server;

#[cfg(test)]
mod tests;

use std::path::PathBuf;

pub use server::{serve, UiServer};

/// Configuration for the local AstroSat web UI server (`astrosat ui`).
#[derive(Debug, Clone)]
pub struct UiConfig {
    pub host: String,
    pub port: u16,
    pub data_dir: PathBuf,
    pub open_browser: bool,
    pub threads: Option<usize>,
}

impl Default for UiConfig {
    fn default() -> Self {
        Self {
            host: "127.0.0.1".to_string(),
            port: 8080,
            data_dir: std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
            open_browser: true,
            threads: None,
        }
    }
}
