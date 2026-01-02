use thiserror::Error;

/// Error type for gxf2bed operations.
#[derive(Debug, Error)]
pub enum Gxf2BedError {
    /// Input extension is missing or not supported.
    #[error("unsupported input extension: {0}")]
    UnsupportedExtension(String),
    /// Failed to build a Rayon thread pool.
    #[error("failed to build thread pool: {0}")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),
    /// Wraps reader errors from genepred.
    #[error("reader error: {0}")]
    Reader(#[from] genepred::reader::ReaderError),
    /// Wraps writer errors from genepred.
    #[error("writer error: {0}")]
    Writer(#[from] genepred::WriterError),
    /// Wraps standard I/O errors.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Result alias for gxf2bed operations.
pub type Result<T> = std::result::Result<T, Gxf2BedError>;
