use std::error::Error;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use colored::Colorize;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use hashbrown::HashMap;
use indoc::indoc;
use rayon::prelude::*;

use crate::cli::Args;
use crate::gxf::{GenePred, GxfRecord, RecordType};

const VERSION: &str = env!("CARGO_PKG_VERSION");

pub fn convert(args: Args) {
    let mut sep = b' ';

    let contents = match args.gxf.extension().and_then(|s| s.to_str()) {
        Some("gz") => {
            match Path::new(args.gxf.file_stem().unwrap())
                .extension()
                .expect("ERROR: No extension found")
                .to_str()
            {
                Some("gff") | Some("gff3") => {
                    sep = b'=';
                }
                _ => (),
            };
            with_gz(&args.gxf).expect("ERROR: Could not read GZ file")
        }
        Some("gtf") => raw(&args.gxf).expect("ERROR: Could not read GTF file"),
        Some("gff") | Some("gff3") => {
            sep = b'=';
            raw(&args.gxf).expect("ERROR: Could not read GFF file")
        }
        _ => panic!("ERROR: Not a GTF/GFF. Wrong file format!"),
    };

    let data = to_bed(&contents, args.parent, args.child, args.feature, sep)
        .expect("ERROR: Could not parse GTF/GFF file");
    log::info!("{} records parsed", data.len());

    write_obj(&args.output, data);
}

pub fn to_bed<'a>(
    content: &str,
    parent: String,
    child: String,
    feature: String,
    sep: u8,
) -> Result<HashMap<String, GenePred>, &'static str> {
    let rs = content
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| match sep {
            b' ' => GxfRecord::parse::<b' '>(row, &feature).ok(),
            b'=' => GxfRecord::parse::<b'='>(row, &feature).ok(),
            _ => None,
        })
        .fold(
            || HashMap::new(),
            |mut acc, record| {
                let feature = record.attr.feature().to_owned();
                let entry = acc.entry(feature).or_insert_with(GenePred::new);

                if record.feature == parent {
                    entry.chr = record.chr.to_owned();
                    entry.start = record.start;
                    entry.end = record.end;
                    entry.strand = record.strand;
                    entry.record_type = RecordType::Parent;
                } else if record.feature == child {
                    entry.chr = record.chr.to_owned();
                    entry.strand = record.strand;
                    entry.start = record.start.min(entry.start);
                    entry.end = record.end.max(entry.end);
                    entry
                        .exons
                        .insert((record.start, record.end - record.start));
                    if entry.record_type != RecordType::Parent {
                        entry.record_type = RecordType::Child;
                    }
                }

                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut left, right| {
                for (feature, info) in right {
                    let entry = left.entry(feature).or_insert_with(GenePred::new);
                    entry.merge(info);
                }
                left
            },
        );

    Ok(rs)
}

pub fn raw<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let mut file = File::open(f)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    Ok(contents)
}

pub fn with_gz<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let file = File::open(f)?;
    let mut decoder = GzDecoder::new(file);
    let mut contents = String::new();

    decoder.read_to_string(&mut contents)?;
    Ok(contents)
}

pub fn max_mem_usage_mb() -> f64 {
    let rusage = unsafe {
        let mut rusage = std::mem::MaybeUninit::uninit();
        libc::getrusage(libc::RUSAGE_SELF, rusage.as_mut_ptr());
        rusage.assume_init()
    };
    let maxrss = rusage.ru_maxrss as f64;
    if cfg!(target_os = "macos") {
        maxrss / 1024.0 / 1024.0
    } else {
        maxrss / 1024.0
    }
}

pub fn write_obj<P: AsRef<Path> + Debug>(filename: P, data: HashMap<String, GenePred>) {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    log::info!("Writing to {:?}", filename);

    let mut writer: Box<dyn Write> = match filename.as_ref().extension() {
        Some(ext) if ext == "gz" => {
            Box::new(BufWriter::new(GzEncoder::new(f, Compression::fast())))
        }
        _ => Box::new(BufWriter::new(f)),
    };

    let mut skips = 0;
    for (transcript, info) in data.into_iter() {
        if info.exons.is_empty() {
            skips += 1;
            continue;
        }

        let (exon_sizes, exon_starts) = info.get_exons_info();
        let (cds_start, cds_end) = info.get_cds();

        if (cds_start >= cds_end) || (info.start >= info.end) {
            log::error!("ERROR: start >= end in record {:?}", info);
            std::process::exit(1);
        }

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            info.chr,
            info.start,
            info.end,
            transcript,
            "0",
            info.strand,
            cds_start,
            cds_end,
            "0",
            info.get_exon_count(),
            exon_sizes,
            exon_starts,
        );
        writeln!(writer, "{}", line).unwrap();
    }

    log::warn!("Skipped {} records with no childs!", skips);
    log::info!("Done writing!");
}

pub fn initialize() {
    println!(
        "{}\n{}\n{}\n",
        "\n##### GXF2BED #####".bright_magenta().bold(),
        indoc!(
            "Fastest GTF/GFF-to-BED converter chilling around.
        Repository: https://github.com/alejandrogzi/gxf2bed
        Feel free to contact the developer if any issue/bug is found."
        ),
        format!("Version: {}", VERSION)
    );
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_to_bed_exon_child() {
        let content = r#"chr1	HAVANA	transcript	92832040	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92832040	92832117	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	start_codon	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92833389	92833458	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92833389	92833458	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	five_prime_utr	92832040	92832114	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	three_prime_utr	92841863	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";"#;

        let data = to_bed(
            &content,
            "transcript".to_string(),
            "exon".to_string(),
            "transcript_id".to_string(),
            b' ',
        )
        .expect("ERROR: Could not parse GTF file");

        assert_eq!(data.len(), 1);
        assert_eq!(data.get("RPL5-202").unwrap().exons.len(), 2);
        assert_eq!(data.get("RPL5-202").unwrap().start, 92832039);
        assert_eq!(data.get("RPL5-202").unwrap().end, 92841924);
        assert_eq!(
            data.get("RPL5-202").unwrap().strand,
            crate::gxf::Strand::Forward
        );
        assert_eq!(
            data.get("RPL5-202").unwrap().record_type,
            crate::gxf::RecordType::Parent
        );
        assert_eq!(data.get("RPL5-202").unwrap().get_exon_count(), 2);
        assert_eq!(
            data.get("RPL5-202").unwrap().get_exons_info(),
            (String::from("78,70,"), String::from("0,1349,"))
        );
    }

    #[test]
    fn test_to_bed_cds_child() {
        let content = r#"chr1	HAVANA	transcript	92832040	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92832040	92832117	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	start_codon	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92833389	92833458	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92833389	92833458	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	five_prime_utr	92832040	92832114	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	three_prime_utr	92841863	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";"#;

        let data = to_bed(
            &content,
            "transcript".to_string(),
            "CDS".to_string(),
            "transcript_id".to_string(),
            b' ',
        )
        .expect("ERROR: Could not parse GTF file");

        assert_eq!(data.len(), 1);
        assert_eq!(data.get("RPL5-202").unwrap().exons.len(), 2);
        assert_eq!(data.get("RPL5-202").unwrap().start, 92832039);
        assert_eq!(data.get("RPL5-202").unwrap().end, 92841924);
        assert_eq!(
            data.get("RPL5-202").unwrap().strand,
            crate::gxf::Strand::Forward
        );
        assert_eq!(
            data.get("RPL5-202").unwrap().record_type,
            crate::gxf::RecordType::Parent
        );
        assert_eq!(data.get("RPL5-202").unwrap().get_exon_count(), 2);
        assert_eq!(
            data.get("RPL5-202").unwrap().get_exons_info(),
            (String::from("3,70,"), String::from("75,1349,"))
        );
    }

    #[test]
    fn test_to_bed_five_utr_child() {
        let content = r#"chr1	HAVANA	transcript	92832040	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92832040	92832117	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	start_codon	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92833389	92833458	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92833389	92833458	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	five_prime_utr	92832040	92832114	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	three_prime_utr	92841863	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";"#;

        let data = to_bed(
            &content,
            "transcript".to_string(),
            "five_prime_utr".to_string(),
            "transcript_id".to_string(),
            b' ',
        )
        .expect("ERROR: Could not parse GTF file");

        assert_eq!(data.len(), 1);
        assert_eq!(data.get("RPL5-202").unwrap().exons.len(), 1);
        assert_eq!(data.get("RPL5-202").unwrap().start, 92832039);
        assert_eq!(data.get("RPL5-202").unwrap().end, 92841924);
        assert_eq!(
            data.get("RPL5-202").unwrap().strand,
            crate::gxf::Strand::Forward
        );
        assert_eq!(
            data.get("RPL5-202").unwrap().record_type,
            crate::gxf::RecordType::Parent
        );
        assert_eq!(data.get("RPL5-202").unwrap().get_exon_count(), 1);
        assert_eq!(
            data.get("RPL5-202").unwrap().get_exons_info(),
            (String::from("75,"), String::from("0,"))
        );
    }

    #[test]
    fn test_to_bed_three_utr_child() {
        let content = r#"chr1	HAVANA	transcript	92832040	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92832040	92832117	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	start_codon	92832115	92832117	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	exon	92833389	92833458	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	CDS	92833389	92833458	.	+	0	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	five_prime_utr	92832040	92832114	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";
        chr1	HAVANA	three_prime_utr	92841863	92841924	.	+	.	gene_symbol "RPL5"; gene_id "ENSG00000122406.14"; gene_name "RPL5"; transcript_id "RPL5-202"; transcript_name "RPL5-202";"#;

        let data = to_bed(
            &content,
            "transcript".to_string(),
            "three_prime_utr".to_string(),
            "transcript_id".to_string(),
            b' ',
        )
        .expect("ERROR: Could not parse GTF file");

        assert_eq!(data.len(), 1);
        assert_eq!(data.get("RPL5-202").unwrap().exons.len(), 1);
        assert_eq!(data.get("RPL5-202").unwrap().start, 92832039);
        assert_eq!(data.get("RPL5-202").unwrap().end, 92841924);
        assert_eq!(
            data.get("RPL5-202").unwrap().strand,
            crate::gxf::Strand::Forward
        );
        assert_eq!(
            data.get("RPL5-202").unwrap().record_type,
            crate::gxf::RecordType::Parent
        );
        assert_eq!(data.get("RPL5-202").unwrap().get_exon_count(), 1);
        assert_eq!(
            data.get("RPL5-202").unwrap().get_exons_info(),
            (String::from("62,"), String::from("9823,"))
        );
    }
}
