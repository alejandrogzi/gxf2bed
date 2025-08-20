//! The fastest GTF/GFF-to-BED converter chilling around
//! Alejandro Gonzales-Irribarren, 2025

use std::collections::BTreeSet;
use std::error::Error;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use colored::Colorize;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use hashbrown::HashMap;
use indoc::indoc;

use crate::cli::Args;
use crate::gxf::{GenePred, GxfRecord, RecordType};

const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Converts a GXF file to a BED file based on the provided arguments.
///
/// This function reads the input GXF file, detects its format (GTF/GFF/GFF3, and gzipped versions),
/// parses the content into `GenePred` records, and then writes them to an output BED file.
/// It uses the `parent`, `child`, and `feature` arguments to define how records are
/// grouped and processed.
///
/// # Arguments
/// * `args` - An `Args` struct containing the input file path, output file path,
///            and parameters for parent, child, and feature identification.
///
/// # Panics
/// * If the input GXF file cannot be read (e.g., file not found, permission denied).
/// * If the file extension is not `gtf`, `gff`, `gff3`, or their gzipped equivalents.
/// * If the `to_bed` function returns an error during parsing.
/// * If `write_obj` encounters an error during file writing.
///
/// # Example
/// ```rust, ignore
/// use std::path::PathBuf;
/// use gxf2bed::{convert, Args}; //
///
/// let args = Args {
///     gxf: PathBuf::from("input.gtf"),
///     output: PathBuf::from("output.bed"),
///     parent: "gene".to_string(),
///     child: "exon".to_string(),
///     feature: "gene_id".to_string(),
/// };
/// convert(args);
/// ```
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

    let data = to_bed(&contents, args.parent, &args.child, args.feature, sep)
        .expect("ERROR: Could not parse GTF/GFF file");
    log::info!("{} records parsed", data.len());

    write_obj(&args.output, data, args.child);
}

/// Parses GXF content (GTF/GFF) into a HashMap of `GenePred` structures.
///
/// This function iterates through lines of GXF content, filters out comments,
/// parses each record using `GxfRecord::parse`, and aggregates them into `GenePred`
/// entries based on a feature key (e.g., gene_id, transcript_id). It differentiates
/// between "parent" and "child" features to construct the `GenePred` structure.
///
/// # Arguments
/// * `content` - A string slice containing the entire GXF file content.
/// * `parent` - The GXF feature type that represents the "parent" entity (e.g., "gene", "transcript").
/// * `child` - The GXF feature type that represents the "child" entity (e.g., "exon", "CDS").
/// * `feature` - The attribute key to use for grouping records into `GenePred` (e.g., "gene_id").
/// * `sep` - The byte separator used within the attributes column (e.g., `b' '` for GTF, `b'='` for GFF).
///
/// # Returns
/// A `Result` containing a `HashMap<String, GenePred>` where keys are the feature IDs
/// and values are the consolidated `GenePred` objects, or a `&'static str` error if parsing fails.
///
/// # Example
/// ```rust, ignore
/// use std::collections::HashMap;
/// use gxf2bed::{to_bed, GenePred, GxfRecord}; //
///
/// let content = "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_id \"geneA\";\n\
///                chr1\t.\texon\t100\t200\t.\t+\t.\tgene_id \"geneA\";\n\
///                chr1\t.\texon\t300\t400\t.\t+\t.\tgene_id \"geneA\";";
///
/// let gene_preds = to_bed(content, "gene".to_string(), &"exon".to_string(), "gene_id".to_string(), b' ').unwrap();
///
/// assert!(gene_preds.contains_key("geneA"));
/// let geneA = gene_preds.get("geneA").unwrap();
/// assert_eq!(geneA.chr, "chr1");
/// assert_eq!(geneA.start, 100 -1); // 0-based
/// assert_eq!(geneA.end, 500);
/// assert_eq!(geneA.exons.len(), 2);
/// ```
pub fn to_bed<'a>(
    content: &str,
    parent: String,
    child: &String,
    feature: String,
    sep: u8,
) -> Result<HashMap<String, GenePred>, &'static str> {
    let mut acc: HashMap<String, GenePred> = HashMap::new();

    content
        .lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|row| match sep {
            b' ' => GxfRecord::parse::<b' '>(row, &feature).ok(),
            b'=' => GxfRecord::parse::<b'='>(row, &feature).ok(),
            _ => None,
        })
        .for_each(|record| {
            let key = record.attr.feature().to_owned();
            let entry = acc.entry(key).or_insert_with(GenePred::new);

            // INFO: parent sets transcript bounds
            if record.feature == parent {
                entry.chr = record.chr.clone();
                entry.start = record.start;
                entry.end = record.end;
                entry.strand = record.strand;
                entry.record_type = RecordType::Parent;
            }

            // INFO: only for CDS is necessary to store exonic blocks!
            if child == "CDS" {
                if record.feature == "exon" {
                    entry
                        .blocks
                        .get_or_insert_with(BTreeSet::new)
                        .insert((record.start, record.end - record.start));
                }
            }

            // INFO: child sets block bounds
            if record.feature == child {
                entry
                    .exons
                    .insert((record.start, record.end - record.start));

                if entry.record_type != RecordType::Parent {
                    entry.record_type = RecordType::Child;
                }
            }
        });

    Ok(acc)
}

/// Reads the entire content of a regular (non-gzipped) file into a String.
///
/// # Type Parameters
/// * `P` - The type of the file path, which must implement `AsRef<Path>` and `Debug`.
///
/// # Arguments
/// * `f` - The path to the file to read.
///
/// # Returns
/// A `Result` containing the file content as a `String` or a `Box<dyn Error>` if an I/O error occurs.
///
/// # Example
/// ```rust, ignore
/// use std::path::Path;
/// use gxf2bed::raw; //
///
/// // Create a dummy file for testing
/// std::fs::write("test_file.txt", "Hello, world!").unwrap();
///
/// let content = raw(Path::new("test_file.txt")).unwrap();
/// assert_eq!(content, "Hello, world!");
///
/// // Clean up
/// std::fs::remove_file("test_file.txt").unwrap();
/// ```
pub fn raw<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let mut file = File::open(f)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    Ok(contents)
}

/// Reads the entire content of a gzipped file into a String.
///
/// This function decompresses the gzipped file and then reads its content.
///
/// # Type Parameters
/// * `P` - The type of the file path, which must implement `AsRef<Path>` and `Debug`.
///
/// # Arguments
/// * `f` - The path to the gzipped file to read.
///
/// # Returns
/// A `Result` containing the decompressed file content as a `String` or a `Box<dyn Error>` if an I/O or decompression error occurs.
///
/// # Example
/// ```rust, ignore
/// use std::path::Path;
/// use gxf2bed::with_gz; //
/// use flate2::{Compression, write::GzEncoder};
/// use std::io::Write;
/// use std::fs::File;
///
/// // Create a dummy gzipped file for testing
/// let mut encoder = GzEncoder::new(File::create("test_file.gz").unwrap(), Compression::default());
/// encoder.write_all(b"Hello, gzipped world!").unwrap();
/// encoder.finish().unwrap();
///
/// let content = with_gz(Path::new("test_file.gz")).unwrap();
/// assert_eq!(content, "Hello, gzipped world!");
///
/// // Clean up
/// std::fs::remove_file("test_file.gz").unwrap();
/// ```
pub fn with_gz<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let file = File::open(f)?;
    let mut decoder = GzDecoder::new(file);
    let mut contents = String::new();

    decoder.read_to_string(&mut contents)?;
    Ok(contents)
}

/// Returns the maximum resident set size (RSS) memory usage of the current process in megabytes.
///
/// This function uses the `getrusage` system call to retrieve resource usage information.
/// The unit of `ru_maxrss` varies by OS (kilobytes on Linux, bytes on macOS), so it adjusts
/// the calculation accordingly.
///
/// # Returns
/// The maximum memory usage in megabytes (`f64`).
///
/// # Safety
/// This function uses `unsafe` block to call the `libc::getrusage` C function.
/// It assumes the `getrusage` call is safe and the `rusage` struct is initialized correctly.
///
/// # Example
/// ```rust, ignore
/// // This example will show memory usage, but its value depends on the runtime environment.
/// use gxf2bed::max_mem_usage_mb; //
///
/// let mem_usage = max_mem_usage_mb();
/// println!("Max memory usage: {:.2} MB", mem_usage);
/// ```
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

/// Writes a `HashMap` of `GenePred` objects to a BED file.
///
/// This function processes the `GenePred` data, formats each entry into a BED line,
/// and writes it to the specified output file. It handles both plain text and gzipped
/// output files based on the file extension. Records with no exons are skipped.
/// Special logic is applied if the `child` feature is "CDS" to adjust coordinates
/// for stop codons and potentially replace exon blocks with `blocks` data.
///
/// # Arguments
/// * `filename` - The path to the output BED file. Can be `.bed` or `.bed.gz`.
/// * `data` - A `HashMap` where keys are transcript/gene IDs and values are `GenePred` objects.
/// * `child` - The name of the child feature (e.g., "exon", "CDS"), which influences
///             how exons are handled, especially for CDS records.
///
/// # Panics
/// * If the output file cannot be created.
/// * If `child` is "CDS" but no `blocks` are found for a transcript, leading to an error exit.
/// * If a record has `start >= end` or `cds_start >= cds_end`, leading to an error exit.
/// * If there's an I/O error during writing to the file.
///
/// # Example
/// ```rust, ignore
/// use std::collections::{HashMap, BTreeSet};
/// use std::path::PathBuf;
/// use gxf2bed::{write_obj, GenePred, Strand, RecordType}; //
///
/// let mut data: HashMap<String, GenePred> = HashMap::new();
/// let mut gp = GenePred::new();
/// gp.chr = "chr1".to_string();
/// gp.start = 100;
/// gp.end = 200;
/// gp.strand = Strand::Forward;
/// gp.exons.insert((100, 150 - 100)); // (start, size)
/// gp.exons.insert((160, 200 - 160));
/// gp.record_type = RecordType::Parent;
/// data.insert("transcript1".to_string(), gp);
///
/// // This will create a file named "output.bed"
/// write_obj(&PathBuf::from("output.bed"), data.clone(), "exon".to_string());
///
/// // Example for CDS processing (will require blocks to be set if "CDS" is child)
/// let mut gp_cds = GenePred::new();
/// gp_cds.chr = "chr1".to_string();
/// gp_cds.start = 1000;
/// gp_cds.end = 2000;
/// gp_cds.strand = Strand::Forward;
/// gp_cds.exons.insert((1000, 100)); // A placeholder exon, will be replaced by blocks
/// let mut cds_blocks = BTreeSet::new();
/// cds_blocks.insert((1000, 50));
/// cds_blocks.insert((1100, 70));
/// gp_cds.blocks = Some(cds_blocks);
/// gp_cds.record_type = RecordType::Parent;
/// data.clear();
/// data.insert("transcript_cds".to_string(), gp_cds);
/// // This will create a file named "output_cds.bed"
/// write_obj(&PathBuf::from("output_cds.bed"), data, "CDS".to_string());
/// ```
pub fn write_obj<P: AsRef<Path> + Debug>(
    filename: P,
    data: HashMap<String, GenePred>,
    child: String,
) {
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
    for (transcript, mut info) in data.into_iter() {
        if info.exons.is_empty() {
            skips += 1;
            continue;
        }

        // INFO: storing original target coords to be used as CDS bounds
        let mut cds_start = info.exons.first().map_or(info.start, |(start, _)| *start);
        let mut cds_end = info
            .exons
            .last()
            .map_or(info.end, |(last_start, size)| last_start + size);

        if child == "CDS" {
            // INFO: replacing exonic blocks with gp.blocks to replicate UTR structure
            if let Some(blocks) = info.blocks.as_mut() {
                info.exons.clear();
                for (start, size) in blocks.iter() {
                    info.exons.insert((*start, *size));
                }

                blocks.clear();
            } else {
                log::error!(
                    "ERROR: No blocks found for transcript {} -> {:?}",
                    transcript,
                    info
                );
                std::process::exit(1);
            }

            // INFO: adjusting CDS start/end to include stop codons
            match info.strand {
                crate::gxf::Strand::Forward => {
                    cds_end += 3;
                }
                crate::gxf::Strand::Reverse => {
                    cds_start -= 3;
                }
                _ => {
                    log::error!("ERROR: Unknown strand in record {:?}", info);
                    std::process::exit(1);
                }
            }
        } else {
            // WARN: neccesary for "exon"?
            // INFO: mutating first/last exon bounds to tx.start and tx.end to make target be represented as CDS
            if let (Some((first_start, first_size)), Some((last_start, last_size))) =
                (info.exons.first().cloned(), info.exons.last().cloned())
            {
                // INFO: remove the old entries
                info.exons.remove(&(first_start, first_size));
                info.exons.remove(&(last_start, last_size));

                // INFO: adjust first and last exons
                let first_end = first_start + first_size;
                let new_first_size = first_end.saturating_sub(info.start);
                let new_last_size = info.end.saturating_sub(last_start);

                info.exons.insert((info.start, new_first_size));
                info.exons.insert((last_start, new_last_size));
            }
        }

        let (exon_sizes, exon_starts) = info.get_exons_info();

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

/// Initializes the application by printing a welcome message and version information to the console.
///
/// This function outputs a stylized banner for "GXF2BED", provides a brief description
/// and repository link, and displays the current version of the tool. It uses the `colored`
/// crate for colorful text and `indoc` for multi-line string formatting.
///
/// # Example
/// ```rust, ignore
/// use gxf2bed::initialize;
/// initialize();
/// ```
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
            &"exon".to_string(),
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
            &"CDS".to_string(),
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
            &"five_prime_utr".to_string(),
            "transcript_id".to_string(),
            b' ',
        )
        .expect("ERROR: Could not parse GTF file");

        dbg!(&data);

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
            &"three_prime_utr".to_string(),
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

    #[test]
    fn test_to_bed_stop_codon_coord() {
        let content = r#"chr7	HAVANA	transcript	43824499	43832030	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	exon	43824499	43824591	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	exon	43825510	43825665	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	exon	43826057	43826117	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	CDS	43826057	43826117	.	+	0	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	start_codon	43826057	43826059	.	+	0	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	exon	43826799	43826955	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	CDS	43826799	43826955	.	+	2	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	exon	43828424	43828671	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	CDS	43828424	43828671	.	+	1	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	exon	43829137	43829273	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	CDS	43829137	43829273	.	+	2	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	exon	43831489	43832030	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	CDS	43831489	43831644	.	+	0	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	stop_codon	43831645	43831647	.	+	0	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	UTR	43824499	43824591	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	UTR	43825510	43825665	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";
            chr7	HAVANA	UTR	43831645	43832030	.	+	.	gene_id \"ENSMUSG00000050063.18\"; transcript_id \"ENSMUST00000107968.9\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"Klk6\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"Klk6-001\"; level 2; protein_id \"ENSMUSP00000103602.3\"; transcript_support_level \"1\"; tag \"NAGNAG_splice_site\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS21184.2\"; havana_gene \"OTTMUSG00000022524.6\"; havana_transcript \"OTTMUST00000053947.4\";"#;

        let data = to_bed(
            &content,
            "transcript".to_string(),
            &"CDS".to_string(),
            "transcript_id".to_string(),
            b' ',
        )
        .expect("ERROR: Could not parse GTF file");

        assert_eq!(data.len(), 1);

        let gp = data.iter().next().unwrap().1;

        assert_eq!(gp.start, 43824498);
        assert_eq!(gp.end, 43832030);
        assert_eq!(gp.strand, crate::gxf::Strand::Forward);
        assert_eq!(gp.record_type, crate::gxf::RecordType::Parent);
    }
}
