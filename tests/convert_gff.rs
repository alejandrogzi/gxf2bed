use gxf2bed::{run, BedType, Config};
use indoc::indoc;
use std::path::{Path, PathBuf};

/// Writes a file to the temporary directory and returns its path.
fn write_temp_file(dir: &Path, name: &str, contents: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, contents).unwrap();
    path
}

/// Converts a small GFF3 to BED12 and validates coordinates.
#[test]
fn convert_gff_to_bed12() {
    let dir = tempfile::tempdir().unwrap();
    let gff = indoc! {"
        chr1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Name=tx1;
        chr1\tsrc\texon\t100\t150\t.\t+\t.\tParent=tx1;
        chr1\tsrc\texon\t180\t200\t.\t+\t.\tParent=tx1;
    "};
    let input_path = write_temp_file(dir.path(), "input.gff3", gff.trim());
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
    assert_eq!(fields.len(), 12);
    assert_eq!(fields[0], "chr1");
    assert_eq!(fields[1], "99");
    assert_eq!(fields[2], "200");
    assert_eq!(fields[3], "tx1");
    assert_eq!(fields[5], "+");
    assert_eq!(fields[9], "2");
    assert_eq!(fields[10], "51,21,");
    assert_eq!(fields[11], "0,80,");
}
