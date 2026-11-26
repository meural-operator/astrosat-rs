//! Error types for the AstroSat CZTI Level-1 -> Level-2 pipeline.

use std::path::PathBuf;
use thiserror::Error;

pub type Result<T> = std::result::Result<T, CztiError>;

#[derive(Debug, Error)]
pub enum CztiError {
    #[error("FITS error: {0}")]
    Fits(#[from] astrosat_fits::FitsError),

    #[error("Core astronomy/timing error: {0}")]
    Core(#[from] astrosat_core::CoreError),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid configuration or parameter: {0}")]
    InvalidConfig(String),

    #[error("Missing required file or CALDB entry: {0}")]
    MissingFile(PathBuf),

    #[error("CALDB lookup failed for {codename} ({detnam}): {reason}")]
    CaldbLookup {
        codename: String,
        detnam: String,
        reason: String,
    },

    #[error("Data decoding error: {0}")]
    Decode(String),

    #[error("Source out of CZTI field of view (|tx|, |ty| <= 5.8 deg): tx={tx:.4} deg, ty={ty:.4} deg")]
    OutOfFov { tx: f64, ty: f64 },

    #[error("No good time intervals remaining after filtering")]
    EmptyGti,

    #[error("Incompatible event file: {0}")]
    IncompatibleInput(String),

    #[error("Processing error in stage {stage}: {message}")]
    Stage { stage: &'static str, message: String },
}
