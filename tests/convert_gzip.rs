use flate2::write::GzEncoder;
use flate2::Compression as GzCompression;
use gxf2bed::{run, BedType, Config};
use indoc::indoc;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Writes gz-compressed contents to a file and returns its path.
fn write_gzip_file(dir: &Path, name: &str, contents: &str) -> PathBuf {
    let mut encoder = GzEncoder::new(Vec::new(), GzCompression::default());
    encoder.write_all(contents.as_bytes()).unwrap();
    let gz = encoder.finish().unwrap();

    let path = dir.join(name);
    std::fs::write(&path, gz).unwrap();
    path
}

/// Converts a gzipped GTF to BED12 using buffered input.
#[test]
fn convert_gzipped_gtf() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = indoc! {"
        chr1\tsrc\ttranscript\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; transcript_name \"tx1\";
        chr1\tsrc\texon\t100\t150\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; exon_number \"1\";
        chr1\tsrc\texon\t180\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; exon_number \"2\";
    "};
    let input_path = write_gzip_file(dir.path(), "input.gtf.gz", gtf.trim());
    let output_path = dir.path().join("output.bed");

    let config = Config {
        input: input_path,
        output: output_path.clone(),
        threads: 2,
        parent_feature: None,
        child_features: None,
        parent_attribute: None,
        child_attribute: None,
        bed_type: BedType::Bed12,
        additional_fields: None,
        chunks: 1024,
    };

    run(&config).unwrap();

    let output = std::fs::read_to_string(&output_path).unwrap();
    let line = output.lines().next().unwrap();
    let fields = line.split('\t').collect::<Vec<_>>();
    assert_eq!(fields[0], "chr1");
    assert_eq!(fields[1], "99");
    assert_eq!(fields[2], "200");
    assert_eq!(fields[3], "tx1");
}
