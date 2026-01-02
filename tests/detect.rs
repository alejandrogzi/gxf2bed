use gxf2bed::detect::{detect_input_kind, Compression, InputFormat};
use std::path::Path;

/// Ensures plain GTF input is detected correctly.
#[test]
fn detect_plain_gtf() {
    let kind = detect_input_kind(Path::new("sample.gtf")).unwrap();
    assert_eq!(kind.format, InputFormat::Gtf);
    assert_eq!(kind.compression, Compression::None);
}

/// Ensures GFF3 input with gzip compression is detected correctly.
#[test]
fn detect_gff3_gz() {
    let kind = detect_input_kind(Path::new("sample.gff3.gz")).unwrap();
    assert_eq!(kind.format, InputFormat::Gff);
    assert_eq!(kind.compression, Compression::Gzip);
}

/// Ensures Zstandard compression is detected correctly.
#[test]
fn detect_gtf_zstd() {
    let kind = detect_input_kind(Path::new("sample.gtf.zst")).unwrap();
    assert_eq!(kind.format, InputFormat::Gtf);
    assert_eq!(kind.compression, Compression::Zstd);
}

/// Rejects unsupported extensions.
#[test]
fn detect_rejects_unknown() {
    assert!(detect_input_kind(Path::new("sample.txt")).is_err());
}
