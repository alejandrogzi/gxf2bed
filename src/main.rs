use gxf2bed::*;

use clap::Parser;
use hashbrown::HashMap;
use log::Level;
use rayon::prelude::*;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Parser, Debug)]
#[clap(
    name = "gxf2bed",
    version = "0.1.0",
    author = "Alejandro Gonzales-Irribarren <jose.gonzalesdezavala1@unmsm.edu.pe>",
    about = "Fastest GTF/GFF-to-BED converter chilling around"
)]
struct Args {
    /// The fastest G{T,F}F-to-BED converter chilling around the world!
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
            return Err(ArgError::InvalidInput(err));
        } else if !self.gxf.extension().unwrap().eq("gff")
            & !self.gxf.extension().unwrap().eq("gtf")
            & !self.gxf.extension().unwrap().eq("gff3")
        {
            let err = format!(
                "file {:?} is not a GTF or GFF3 file, please specify the correct format",
                self.gxf
            );
            return Err(ArgError::InvalidInput(err));
        } else if std::fs::metadata(&self.gxf).unwrap().len() == 0 {
            let err = format!("file {:?} is empty", self.gxf);
            return Err(ArgError::InvalidInput(err));
        } else {
            Ok(())
        }
    }

    /// Checks the output file for validity. If the file is not a BED file, an error is returned.
    fn check_output(&self) -> Result<(), ArgError> {
        if !self.output.extension().unwrap().eq("bed") {
            let err = format!("file {:?} is not a BED file", self.output);
            return Err(ArgError::InvalidOutput(err));
        } else {
            Ok(())
        }
    }

    /// Checks the number of threads for validity. The number of threads must be greater than 0
    /// and less than or equal to the number of logical CPUs.
    fn check_threads(&self) -> Result<(), ArgError> {
        if self.threads == 0 {
            let err = format!("number of threads must be greater than 0");
            return Err(ArgError::InvalidThreads(err));
        } else if self.threads > num_cpus::get() {
            let err = format!(
                "number of threads must be less than or equal to the number of logical CPUs"
            );
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

    let contents = reader(&args.gxf)?;
    let data = to_bed(&contents, args.parent, args.child, args.feature)?;
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
    s: &'a str,
    parent: String,
    child: String,
    feature: String,
) -> Result<HashMap<String, HashMap<&str, String>>, &'static str> {
    s.par_lines()
        .map(|line| {
            if !line.starts_with("#") {
                Some(GxfRecord::parse(line, &feature))
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

            if record.feat == parent {
                // args.parent
                entry.insert("chr", record.chr.to_owned());
                entry.insert("start", record.start.to_string());
                entry.insert("end", record.end.to_string());
                entry.insert("strand", record.strand.to_string());
            } else if record.feat == child {
                // args.child
                entry.entry("exons").or_default().push('.');

                let exon_starts = entry.entry("exon_starts").or_insert(String::from(""));
                exon_starts.push_str(&record.start.to_string());
                exon_starts.push_str(",");

                let exon_sizes = entry.entry("exon_sizes").or_insert(String::from(""));
                exon_sizes.push_str(&(record.end - record.start).to_string());
                exon_sizes.push_str(",");
            } else if record.feat == "start_codon" {
                entry.insert("start_codon", record.start.to_string());
            } else if record.feat == "stop_codon" {
                entry.insert("stop_codon", record.start.to_string());
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
