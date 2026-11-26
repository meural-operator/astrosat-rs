use thiserror::Error;

/// Error type for all FITS reading, writing, and validation operations.
#[derive(Debug, Error)]
pub enum FitsError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("invalid FITS header: {0}")]
    InvalidHeader(String),

    #[error("invalid FITS block alignment: {0}")]
    InvalidBlockAlignment(String),

    #[error("unsupported BITPIX value: {0}")]
    UnsupportedBitpix(i64),

    #[error("invalid TFORM specification '{0}'")]
    InvalidTform(String),

    #[error("HDU not found: {0}")]
    HduNotFound(String),

    #[error("column not found: '{0}'")]
    ColumnNotFound(String),

    #[error("column or HDU type mismatch: {0}")]
    TypeMismatch(String),

    #[error("index out of bounds: {0}")]
    OutOfBounds(String),

    #[error("checksum mismatch in HDU {hdu}: expected {expected}, got {actual}")]
    ChecksumMismatch {
        hdu: usize,
        expected: String,
        actual: String,
    },

    #[error("FITS validation error: {0}")]
    ValidationError(String),

    #[error("FITS format error: {0}")]
    FormatError(String),
}

pub type Result<T> = std::result::Result<T, FitsError>;
