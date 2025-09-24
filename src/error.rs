//! Error types for the FastGA-RS library.

use std::path::PathBuf;
use thiserror::Error;

/// Result type alias for FastGA operations.
pub type Result<T> = std::result::Result<T, FastGAError>;

/// Errors that can occur during FastGA operations.
#[derive(Error, Debug)]
pub enum FastGAError {
    /// Input file not found
    #[error("File not found: {0}")]
    FileNotFound(PathBuf),

    /// I/O error during file operations
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),

    /// FastGA execution failed
    #[error("FastGA execution failed: {0}")]
    FastGAExecutionFailed(String),

    /// Failed to parse PAF output
    #[error("Failed to parse PAF output: {0}")]
    PafParseError(String),

    /// Failed to parse CIGAR string
    #[error("Failed to parse CIGAR string: {0}")]
    CigarParseError(String),

    /// Invalid configuration parameter
    #[error("Invalid configuration: {0}")]
    InvalidConfig(String),

    /// FFI error when calling C code
    #[error("FFI error: {0}")]
    FfiError(String),

    /// Alignment operation was cancelled
    #[error("Alignment operation was cancelled")]
    Cancelled,

    /// Memory allocation failed
    #[error("Memory allocation failed: {0}")]
    AllocationError(String),

    /// Invalid sequence format
    #[error("Invalid sequence format: {0}")]
    InvalidSequence(String),

    /// Temporary directory creation failed
    #[error("Failed to create temporary directory")]
    TempDirError,

    /// UTF-8 conversion error
    #[error("UTF-8 conversion error: {0}")]
    Utf8Error(#[from] std::string::FromUtf8Error),

    /// Generic error with custom message
    #[error("{0}")]
    Other(String),
}

impl From<tempfile::PersistError> for FastGAError {
    fn from(_: tempfile::PersistError) -> Self {
        FastGAError::TempDirError
    }
}