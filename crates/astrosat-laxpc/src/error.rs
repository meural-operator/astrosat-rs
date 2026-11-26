//! Error types for the AstroSat LAXPC Level-1 to Level-2 pipeline.

use std::path::PathBuf;
use thiserror::Error;

/// Unified error type for `astrosat-laxpc`.
#[derive(Debug, Error)]
pub enum LaxpcError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FITS error: {0}")]
    Fits(String),

    #[error("Core science/timing error: {0}")]
    Core(String),

    #[error("Invalid LAXPC unit identifier: {0} (expected 1/lxp10, 2/lxp20, or 3/lxp30)")]
    InvalidUnit(String),

    #[error("Invalid configuration: {0}")]
    Config(String),

    #[error("Frame parsing/sorting error: {0}")]
    Frame(String),

    #[error("TCT conversion error: {0}")]
    Tct(String),

    #[error("GTI generation error: {0}")]
    Gti(String),

    #[error("Background/gain model error: {0}")]
    Background(String),

    #[error("Channel map error in {path}: {message}")]
    ChannelMap { path: PathBuf, message: String },

    #[error("Missing required input file: {0}")]
    MissingInput(PathBuf),
}

/// Convenience result alias for LAXPC operations.
pub type Result<T> = std::result::Result<T, LaxpcError>;
