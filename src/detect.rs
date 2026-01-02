use crate::error::{Gxf2BedError, Result};
use std::path::{Path, PathBuf};

/// Supported input formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    Gtf,
    Gff,
}

/// Supported compression formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    None,
    Gzip,
    Zstd,
    Bzip2,
}

impl Compression {
    /// Returns true when the input is compressed.
    /// 
    /// # Returns
    /// 
    /// Returns true if the compression type is not `Compression::None`.
    /// 
    /// # Example
    /// 
    /// ```rust, ignore
    /// use gxf2bed::detect::Compression;
    /// 
    /// assert!(Compression::Gzip.is_compressed());
    /// assert!(!Compression::None.is_compressed());
    /// ```
    pub fn is_compressed(self) -> bool {
        !matches!(self, Compression::None)
    }
}

/// Describes the detected input kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InputKind {
    /// Input format (GTF or GFF).
    pub format: InputFormat,
    /// Compression format.
    pub compression: Compression,
}

/// Detects input format and compression from the file extension(s).
/// 
/// Analyzes file extensions to determine both the format (GTF/GFF)
/// and compression type. Handles nested extensions like `.gtf.gz`.
/// 
/// # Arguments
/// 
/// * `path` - Path to the input file
/// 
/// # Returns
/// 
/// Returns InputKind containing detected format and compression.
/// 
/// # Errors
/// 
/// Returns an error if the file extension is not supported.
/// 
/// # Example
/// 
/// ```rust, ignore
/// use gxf2bed::detect::detect_input_kind;
/// use std::path::Path;
/// 
/// let kind = detect_input_kind(Path::new("file.gtf.gz"))?;
/// // Returns InputKind with Gtf format and Gzip compression
/// ```
pub fn detect_input_kind(path: &Path) -> Result<InputKind> {
    let ext = extension_lowercase(path)
        .ok_or_else(|| Gxf2BedError::UnsupportedExtension(path.display().to_string()))?;

    if let Some(compression) = compression_from_extension(&ext) {
        let inner_ext = nested_extension(path)
            .ok_or_else(|| Gxf2BedError::UnsupportedExtension(path.display().to_string()))?;
        let format = format_from_extension(&inner_ext)
            .ok_or_else(|| Gxf2BedError::UnsupportedExtension(path.display().to_string()))?;
        return Ok(InputKind {
            format,
            compression,
        });
    }

    let format = format_from_extension(&ext)
        .ok_or_else(|| Gxf2BedError::UnsupportedExtension(path.display().to_string()))?;
    Ok(InputKind {
        format,
        compression: Compression::None,
    })
}

/// Extracts the lowercase extension from a path.
/// 
/// Gets the file extension and converts it to lowercase for
/// case-insensitive comparison.
/// 
/// # Arguments
/// 
/// * `path` - Path to extract extension from
/// 
/// # Returns
/// 
/// Returns the lowercase extension as a String, or None if no extension exists.
/// 
/// # Example
/// 
/// ```rust, ignore
/// use gxf2bed::detect::extension_lowercase;
/// use std::path::Path;
/// 
/// let ext = extension_lowercase(Path::new("file.GTF"));
/// assert_eq!(ext, Some("gtf".to_string()));
/// ```
fn extension_lowercase(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.to_ascii_lowercase())
}

/// Determines the input format based on extension.
/// 
/// Maps file extensions to their corresponding InputFormat variants.
/// Supports .gtf, .gff, and .gff3 extensions.
/// 
/// # Arguments
/// 
/// * `ext` - File extension string (without the dot)
/// 
/// # Returns
/// 
/// Returns the corresponding InputFormat, or None if extension is unsupported.
/// 
/// # Example
/// 
/// ```rust, ignore
/// use gxf2bed::detect::format_from_extension;
/// use gxf2bed::detect::InputFormat;
/// 
/// let format = format_from_extension("gtf");
/// assert_eq!(format, Some(InputFormat::Gtf));
/// ```
fn format_from_extension(ext: &str) -> Option<InputFormat> {
    match ext {
        "gtf" => Some(InputFormat::Gtf),
        "gff" | "gff3" => Some(InputFormat::Gff),
        _ => None,
    }
}

/// Determines compression based on extension.
/// 
/// Maps file extensions to their corresponding Compression variants.
/// Supports common compression formats like gzip, zstd, and bzip2.
/// 
/// # Arguments
/// 
/// * `ext` - File extension string (without the dot)
/// 
/// # Returns
/// 
/// Returns the corresponding Compression variant, or None if not a compression format.
/// 
/// # Example
/// 
/// ```rust, ignore
/// use gxf2bed::detect::compression_from_extension;
/// use gxf2bed::detect::Compression;
/// 
/// let compression = compression_from_extension("gz");
/// assert_eq!(compression, Some(Compression::Gzip));
/// ```
fn compression_from_extension(ext: &str) -> Option<Compression> {
    match ext {
        "gz" | "gzip" => Some(Compression::Gzip),
        "zst" | "zstd" => Some(Compression::Zstd),
        "bz2" | "bzip2" => Some(Compression::Bzip2),
        _ => None,
    }
}

/// Returns the inner extension for compressed files (e.g., `.gtf.gz` -> `gtf`).
/// 
/// For files with nested extensions like compressed GTF/GFF files,
/// extracts the inner file format extension.
/// 
/// # Arguments
/// 
/// * `path` - Path to the compressed file
/// 
/// # Returns
/// 
/// Returns the inner extension (format type), or None if not found.
/// 
/// # Example
/// 
/// ```rust, ignore
/// use gxf2bed::detect::nested_extension;
/// use std::path::Path;
/// 
/// let inner_ext = nested_extension(Path::new("file.gtf.gz"));
/// assert_eq!(inner_ext, Some("gtf".to_string()));
/// ```
fn nested_extension(path: &Path) -> Option<String> {
    let stem = path.file_stem()?.to_str()?;
    extension_lowercase(&PathBuf::from(stem))
}
