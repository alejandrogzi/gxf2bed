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

/// Converts GFF3 using a gene parent and validates that ID wins over Name.
#[test]
fn convert_gff_gene_parent_uses_configured_parent_attribute() {
    let dir = tempfile::tempdir().unwrap();
    let gff = indoc! {"
        scf_000005_A_fumigatus_A1163\tVEuPathDB\tgene\t1488183\t1490405\t.\t+\t.\tID=AFUB_068240;Name=prm1;description=Has domain(s) with predicted role in plasma membrane fusion involved in cytogamy and mating projection tip localization
        scf_000005_A_fumigatus_A1163\tVEuPathDB\tmRNA\t1488183\t1490405\t.\t+\t.\tID=AFUB_068240-T;Parent=AFUB_068240;description=Has domain(s) with predicted role in plasma membrane fusion involved in cytogamy and mating projection tip localization
        scf_000003_A_fumigatus_A1163\tVEuPathDB\tgene\t1802685\t1806924\t.\t+\t.\tID=AFUB_040320;description=Has domain(s) with predicted binding activity
        scf_000003_A_fumigatus_A1163\tVEuPathDB\tmRNA\t1802685\t1806924\t.\t+\t.\tID=AFUB_040320-T;Parent=AFUB_040320;description=Has domain(s) with predicted binding activity
    "};
    let input_path = write_temp_file(dir.path(), "input.gff3", gff.trim());
    let output_path = dir.path().join("output.bed");

    let config = Config {
        input: input_path,
        output: output_path.clone(),
        threads: 2,
        parent_feature: Some("gene".to_string()),
        child_features: Some(vec!["mRNA".to_string()]),
        parent_attribute: Some("ID".to_string()),
        child_attribute: Some("Parent".to_string()),
        bed_type: BedType::Bed12,
        additional_fields: None,
        chunks: 1024,
    };

    run(&config).unwrap();

    let output = std::fs::read_to_string(&output_path).unwrap();
    let mut names = output
        .lines()
        .map(|line| line.split('\t').nth(3).unwrap())
        .collect::<Vec<_>>();
    names.sort_unstable();

    assert_eq!(names, vec!["AFUB_040320", "AFUB_068240"]);
}

/// Converts GFF3 using an mRNA parent and validates that ID wins over Name.
#[test]
fn convert_gff_mrna_parent_uses_configured_parent_attribute() {
    let dir = tempfile::tempdir().unwrap();
    let gff = indoc! {"
        chr1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Name=display_name;
        chr1\tsrc\texon\t100\t150\t.\t+\t.\tParent=tx1;
        chr1\tsrc\texon\t180\t200\t.\t+\t.\tParent=tx1;
    "};
    let input_path = write_temp_file(dir.path(), "input.gff3", gff.trim());
    let output_path = dir.path().join("output.bed");

    let config = Config {
        input: input_path,
        output: output_path.clone(),
        threads: 2,
        parent_feature: Some("mRNA".to_string()),
        child_features: Some(vec!["exon".to_string()]),
        parent_attribute: Some("ID".to_string()),
        child_attribute: Some("Parent".to_string()),
        bed_type: BedType::Bed12,
        additional_fields: None,
        chunks: 1024,
    };

    run(&config).unwrap();

    let output = std::fs::read_to_string(&output_path).unwrap();
    let name = output.lines().next().unwrap().split('\t').nth(3).unwrap();

    assert_eq!(name, "tx1");
}
