use gxf2bed::*;

use clap::Parser;
use hashbrown::HashMap;
use log::Level;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Parser, Debug)]
#[clap(
    name = "gxf2bed",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "fastest GTF/GFF-to-BED converter chilling around"
)]
struct Args {
    /// fastest G{T,F}F-to-BED converter chilling around the world!
    ///
    /// This program converts GTF/GFF3 files to BED format blazingly fast.
    /// Start by providing the path to the GTF/GFF3 file with -i/--input file.gtf
    /// or -i/--input file.gff3.
    #[clap(
        short = 'i',
        long = "input",
        help = "Path to GTF/GFF file",
        value_name = "GXF",
        required = true
    )]
    gxf: PathBuf,

    /// Output filepath; non-required argument.
    ///
    /// The output file will be a BED file with the same name as the input file.
    #[clap(
        short = 'o',
        long = "output",
        help = "Path to output BED file",
        value_name = "BED",
        required = true
    )]
    output: PathBuf,

    /// Number of threads to use; default is the number of logical CPUs.
    #[clap(
        short = 't',
        long,
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    threads: usize,

    /// Parent feature; default is "transcript_id".
    #[clap(
        short = 'p',
        long = "parent",
        help = "Parent feature",
        value_name = "PARENT",
        default_value = "transcript"
    )]
    parent: String,

    /// Child feature; default is "exon".
    #[clap(
        short = 'c',
        long = "child",
        help = "Child feature",
        value_name = "CHILD",
        default_value = "exon"
    )]
    child: String,

    /// Feature to extract; default is "transcript_id".
    #[clap(
        short = 'f',
        long = "feature",
        help = "Feature to extract",
        value_name = "FEATURE",
        default_value = "transcript_id"
    )]
    feature: String,
}

impl Args {
    /// Checks all the arguments for validity using validate_args()
    pub fn check(&self) -> Result<(), ArgError> {
        self.validate_args()
    }

    /// Checks the input file for validity. The file must exist and be a GTF or GFF3 file.
    /// If the file does not exist, an error is returned.
    fn check_input(&self) -> Result<(), ArgError> {
        if !self.gxf.exists() {
            let err = format!("file {:?} does not exist", self.gxf);
            Err(ArgError::InvalidInput(err))
        } else if std::fs::metadata(&self.gxf).unwrap().len() == 0 {
            let err = format!("file {:?} is empty", self.gxf);
            return Err(ArgError::InvalidInput(err));
        } else {
            Ok(())
        }
    }

    /// Checks the output file for validity. If the file is not a BED file, an error is returned.
    fn check_output(&self) -> Result<(), ArgError> {
        if !self.output.extension().unwrap().eq("bed") & !self.output.extension().unwrap().eq("gz")
        {
            let err = format!("file {:?} is not a BED file", self.output);
            Err(ArgError::InvalidOutput(err))
        } else {
            Ok(())
        }
    }

    /// Checks the number of threads for validity. The number of threads must be greater than 0
    /// and less than or equal to the number of logical CPUs.
    fn check_threads(&self) -> Result<(), ArgError> {
        if self.threads == 0 {
            let err = "number of threads must be greater than 0".to_string();
            Err(ArgError::InvalidThreads(err))
        } else if self.threads > num_cpus::get() {
            let err = "number of threads must be less than or equal to the number of logical CPUs"
                .to_string();
            return Err(ArgError::InvalidThreads(err));
        } else {
            Ok(())
        }
    }

    /// Validates all the arguments
    fn validate_args(&self) -> Result<(), ArgError> {
        self.check_input()?;
        self.check_output()?;
        self.check_threads()?;
        Ok(())
    }
}

#[derive(Debug, Error)]
pub enum ArgError {
    /// The input file does not exist or is not a GTF or GFF3 file.
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    /// The output file is not a BED file.
    #[error("Invalid output: {0}")]
    InvalidOutput(String),

    /// The number of threads is invalid.
    #[error("Invalid number of threads: {0}")]
    InvalidThreads(String),
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    msg();
    simple_logger::init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();
    args.check().unwrap_or_else(|e| {
        log::error!("{}", e);
        std::process::exit(1);
    });
    log::info!("{:?}", args);

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    log::info!("Using {} threads", args.threads);

    let st = std::time::Instant::now();
    let start_mem = max_mem_usage_mb();
    let mut sep = b' ';

    let contents = match args.gxf.extension().and_then(|s| s.to_str()) {
        Some("gz") => {
            match Path::new(args.gxf.file_stem().unwrap())
                .extension()
                .expect("ERROR: No extension found")
                .to_str()
            {
                Some("gff") | Some("gff3") => {
                    sep = b';';
                }
                _ => (),
            };
            with_gz(&args.gxf)?
        }
        Some("gtf") => raw(&args.gxf)?,
        Some("gff") | Some("gff3") => {
            sep = b'=';
            raw(&args.gxf)?
        }
        _ => panic!("ERROR: Not a GTF/GFF. Wrong file format!"),
    };

    let data = to_bed(&contents, args.parent, args.child, args.feature, sep)?;

    log::info!("{} records parsed", data.len());

    let mut liner = data
        .into_par_iter()
        .filter(|(_, v)| !v.is_empty())
        .collect::<Vec<_>>();

    liner.par_sort_unstable_by(|a, b| {
        let a_chr = a.1.get("chr").unwrap().parse::<String>().unwrap();
        let b_chr = b.1.get("chr").unwrap().parse::<String>().unwrap();
        let a_start = a.1["start"].parse::<i32>().unwrap();
        let b_start = b.1["start"].parse::<i32>().unwrap();
        a_chr.cmp(&b_chr).then(a_start.cmp(&b_start))
    });

    write_obj(&args.output, liner);

    let elapsed = st.elapsed();
    let mem = (max_mem_usage_mb() - start_mem).max(0.0);

    log::info!("Elapsed: {:.4?} secs", elapsed.as_secs_f32());
    log::info!("Memory: {:.2} MB", mem);

    Ok(())
}

pub fn to_bed<'a>(
    s: &str,
    parent: String,
    child: String,
    feature: String,
    sep: u8,
) -> Result<HashMap<String, HashMap<&'a str, String>>, &'static str> {
    s.par_lines()
        .map(|line| {
            if !line.starts_with("#") {
                match sep {
                    b' ' => Some(GxfRecord::parse::<b' '>(line, &feature)),
                    b'=' => Some(GxfRecord::parse::<b'='>(line, &feature)),
                    _ => None,
                }
            } else {
                None
            }
        })
        .filter_map(|x| x)
        .try_fold_with(HashMap::new(), |mut acc, record| {
            let record = record?;

            let tx_id = if !record.attr.feature().is_empty() {
                record.attr.feature().to_owned()
            } else {
                // continue with the next record
                return Ok(acc);
            };

            let entry = acc.entry(tx_id).or_insert(HashMap::new());

            if !parent.is_empty() {
                send_ft(&record.feat, entry, &record, &parent, &child);
            } else {
                entry.insert("chr", record.chr.to_owned());
                entry.insert("start", record.start.to_string());
                entry.insert("end", record.end.to_string());
                entry.insert("strand", record.strand.to_string());
                to_exon(entry, &record);
            }

            Ok(acc)
        }) // end fold
        .try_reduce_with(|mut map1, map2| {
            for (k, v) in map2 {
                let entry = map1.entry(k).or_insert(HashMap::new());
                for (k2, v2) in v {
                    entry.insert(k2, v2);
                }
            }
            Ok(map1)
        }) // end reduce
        .unwrap_or(Err("Error converting GTF/GFF3 to BED"))
}

fn send_ft(
    feat: &str,
    entry: &mut HashMap<&str, String>,
    record: &GxfRecord,
    parent: &str,
    child: &str,
) {
    if feat == parent {
        entry.insert("chr", record.chr.to_owned());
        entry.insert("start", record.start.to_string());
        entry.insert("end", record.end.to_string());
        entry.insert("strand", record.strand.to_string());
    } else if feat == child {
        to_exon(entry, record);
    } else {
        match feat {
            "start_codon" => {
                entry.insert("start_codon", record.start.to_string());
            }
            "stop_codon" => {
                entry.insert("stop_codon", record.start.to_string());
            }
            _ => {}
        }
    }
}

fn to_exon(entry: &mut HashMap<&str, String>, record: &GxfRecord) {
    entry.entry("exons").or_default().push('.');

    let exon_starts = entry.entry("exon_starts").or_insert(String::from(""));
    exon_starts.push_str(&record.start.to_string());
    exon_starts.push(',');

    let exon_sizes = entry.entry("exon_sizes").or_insert(String::from(""));
    exon_sizes.push_str(&(record.end - record.start).to_string());
    exon_sizes.push(',');
}
