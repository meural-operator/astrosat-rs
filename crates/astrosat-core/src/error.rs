use thiserror::Error;

/// Error type for `astrosat-core` timing, TCT, GTI, MKF, orbit, barycentric, and CALDB operations.
#[derive(Debug, Error)]
pub enum CoreError {
    #[error("FITS error: {0}")]
    Fits(#[from] astrosat_fits::FitsError),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("time {0} out of valid range")]
    TimeOutOfRange(f64),

    #[error("invalid TCT table: {0}")]
    InvalidTct(String),

    #[error("invalid orbit table: {0}")]
    InvalidOrbit(String),

    #[error("invalid MKF table: {0}")]
    InvalidMkf(String),

    #[error("invalid ephemeris: {0}")]
    InvalidEphemeris(String),

    #[error("CALDB error: {0}")]
    CaldbError(String),

    #[error("missing required keyword: {0}")]
    MissingKeyword(String),

    #[error("invalid input: {0}")]
    InvalidInput(String),
}

pub type Result<T> = std::result::Result<T, CoreError>;
