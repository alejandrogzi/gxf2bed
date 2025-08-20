//! The fastest GTF/GFF-to-BED converter chilling around
//! Alejandro Gonzales-Irribarren, 2025

use clap::Parser;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Parser, Debug)]
#[clap(
    name = "gxf2bed",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "fastest GTF/GFF-to-BED converter chilling around"
)]
pub struct Args {
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
    pub gxf: PathBuf,

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
    pub output: PathBuf,

    /// Number of threads to use; default is the number of logical CPUs.
    #[clap(
        short = 't',
        long,
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    /// Parent feature; default is "transcript_id".
    #[clap(
        short = 'p',
        long = "parent",
        help = "Parent feature",
        value_name = "PARENT",
        default_value = "transcript"
    )]
    pub parent: String,

    /// Child feature; default is "exon".
    #[clap(
        short = 'c',
        long = "child",
        help = "Child feature",
        value_name = "CHILD",
        default_value = "exon"
    )]
    pub child: String,

    /// Feature to extract; default is "transcript_id".
    #[clap(
        short = 'f',
        long = "feature",
        help = "Feature to extract",
        value_name = "FEATURE",
        default_value = "transcript_id"
    )]
    pub feature: String,
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
