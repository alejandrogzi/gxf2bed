use clap::Parser;
use std::path::PathBuf;

/// CLI arguments for gxf2bed.
#[derive(Parser, Debug)]
#[clap(
    name = "gxf2bed",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "Fastest GTF/GFF-to-BED converter chilling around"
)]
pub struct Args {
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
    pub gxf: PathBuf,

    /// Output filepath; required argument.
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
        short = 'T',
        long,
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    /// Parent feature; default is "transcript_id".
    #[clap(
        short = 'F',
        long = "parent-feature",
        help = "Parent feature",
        value_name = "PARENT"
    )]
    pub parent_feature: Option<String>,

    /// Child feature; default is "exon".
    #[clap(
        short = 'f',
        long = "child-features",
        help = "Child features",
        value_name = "CHILDS",
        value_delimiter = ',',
        num_args = 1..,
    )]
    pub child_features: Option<Vec<String>>,

    /// Feature to extract.
    #[clap(
        short = 'A',
        long = "parent-attribute",
        help = "Feature to extract",
        value_name = "FEATURE"
    )]
    pub parent_attribute: Option<String>,

    /// Child feature to extract.
    #[clap(
        short = 'a',
        long = "child-attribute",
        help = "Child feature to extract",
        value_name = "CHILD"
    )]
    pub child_attribute: Option<String>,

    /// BED type format.
    #[clap(
        short = 't',
        long = "type",
        help = "BED type format",
        value_name = "BED_TYPE",
        default_value_t = BedType::Bed12
    )]
    pub bed_type: BedType,

    /// BED additional fields (will use GTF/GFF tags).
    #[clap(
        short = 'd',
        long = "additional-fields",
        help = "BED additional fields",
        value_name = "ADDITIONAL",
        value_delimiter = ',',
        num_args = 1..,
    )]
    pub additional_fields: Option<Vec<String>>,

    /// Chunk size for parallel processing.
    #[clap(
        short = 'c',
        long = "chunks",
        help = "Chunk size for parallel processing",
        value_name = "CHUNKS",
        default_value_t = 15000
    )]
    pub chunks: usize,
}

/// Supported output BED formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BedType {
    Bed3,
    Bed4,
    Bed5,
    Bed6,
    Bed9,
    Bed12,
}

impl Default for BedType {
    /// Returns the default BED output format (BED12).
    ///
    /// # Returns
    ///
    /// The default BED12 format variant.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::BedType;
    /// let default = BedType::default();
    /// assert_eq!(default, BedType::Bed12);
    /// ```
    fn default() -> Self {
        BedType::Bed12
    }
}

impl std::str::FromStr for BedType {
    type Err = String;

    /// Parses a BED type from its numeric representation.
    ///
    /// # Arguments
    ///
    /// * `s` - String slice containing the numeric BED type (3, 4, 5, 6, 9, or 12)
    ///
    /// # Returns
    ///
    /// Returns the corresponding BedType variant or an error for invalid values.
    ///
    /// # Errors
    ///
    /// Returns an error if the input string doesn't match a supported BED type.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::BedType;
    /// use std::str::FromStr;
    ///
    /// let bed_type = BedType::from_str("12").unwrap();
    /// assert_eq!(bed_type, BedType::Bed12);
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "3" => Ok(BedType::Bed3),
            "4" => Ok(BedType::Bed4),
            "5" => Ok(BedType::Bed5),
            "6" => Ok(BedType::Bed6),
            "9" => Ok(BedType::Bed9),
            "12" => Ok(BedType::Bed12),
            _ => Err(format!("Invalid BED type: {}", s)),
        }
    }
}

impl std::fmt::Display for BedType {
    /// Formats the BED type as its numeric representation.
    ///
    /// # Arguments
    ///
    /// * `f` - Formatter to write the BED type to
    ///
    /// # Returns
    ///
    /// Returns a formatting result indicating success or failure.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::BedType;
    ///
    /// let bed_type = BedType::Bed12;
    /// assert_eq!(format!("{}", bed_type), "12");
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BedType::Bed3 => write!(f, "3"),
            BedType::Bed4 => write!(f, "4"),
            BedType::Bed5 => write!(f, "5"),
            BedType::Bed6 => write!(f, "6"),
            BedType::Bed9 => write!(f, "9"),
            BedType::Bed12 => write!(f, "12"),
        }
    }
}
