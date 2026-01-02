use gxf2bed::{run, BedType, Config};
use indoc::indoc;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Writes a file to the temporary directory and returns its path.
fn write_temp_file(dir: &Path, name: &str, contents: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, contents).unwrap();
    path
}

/// Converts a small GTF to BED12 and validates coordinates and blocks.
#[test]
fn convert_gtf_to_bed12() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = indoc! {"
        chr1\tsrc\ttranscript\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; transcript_name \"tx1\";
        chr1\tsrc\texon\t100\t150\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; exon_number \"1\";
        chr1\tsrc\texon\t180\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; exon_number \"2\";
        chr2\tsrc\ttranscript\t1000\t1100\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\"; transcript_name \"tx2\";
        chr2\tsrc\texon\t1000\t1050\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\"; exon_number \"1\";
        chr2\tsrc\texon\t1070\t1100\t.\t-\t.\tgene_id \"g2\"; transcript_id \"tx2\"; exon_number \"2\";
    "};
    let input_path = write_temp_file(dir.path(), "input.gtf", gtf.trim());
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
    let mut by_name: HashMap<String, Vec<String>> = HashMap::new();
    for line in output.lines() {
        let fields = line.split('\t').map(|s| s.to_string()).collect::<Vec<_>>();
        assert_eq!(fields.len(), 12);
        by_name.insert(fields[3].clone(), fields);
    }

    let tx1 = by_name.get("tx1").unwrap();
    assert_eq!(tx1[0], "chr1");
    assert_eq!(tx1[1], "99");
    assert_eq!(tx1[2], "200");
    assert_eq!(tx1[5], "+");
    assert_eq!(tx1[9], "2");
    assert_eq!(tx1[10], "51,21,");
    assert_eq!(tx1[11], "0,80,");

    let tx2 = by_name.get("tx2").unwrap();
    assert_eq!(tx2[0], "chr2");
    assert_eq!(tx2[1], "999");
    assert_eq!(tx2[2], "1100");
    assert_eq!(tx2[5], "-");
    assert_eq!(tx2[9], "2");
    assert_eq!(tx2[10], "51,31,");
    assert_eq!(tx2[11], "0,70,");
}
