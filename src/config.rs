use crate::cli::{Args, BedType};
use std::path::PathBuf;

/// Normalized configuration for a conversion run.
#[derive(Clone, Debug)]
pub struct Config {
    /// Input GTF/GFF path.
    pub input: PathBuf,
    /// Output BED path.
    pub output: PathBuf,
    /// Number of threads to use.
    pub threads: usize,
    /// Parent feature override.
    pub parent_feature: Option<String>,
    /// Child feature overrides.
    pub child_features: Option<Vec<String>>,
    /// Parent attribute override.
    pub parent_attribute: Option<String>,
    /// Child attribute override.
    pub child_attribute: Option<String>,
    /// Output BED type.
    pub bed_type: BedType,
    /// Additional fields to include in output.
    pub additional_fields: Option<Vec<String>>,
    /// Chunk size for parallel processing.
    pub chunks: usize,
}

impl Config {
    /// Builds a conversion config from CLI arguments.
    /// 
    /// # Arguments
    /// 
    /// * `args` - Command-line arguments to convert into a configuration
    /// 
    /// # Returns
    /// 
    /// Returns a new Config instance with values copied from the CLI arguments.
    /// 
    /// # Example
    /// 
    /// ```rust, ignore
    /// use gxf2bed::{Args, Config};
    /// use std::path::PathBuf;
    /// 
    /// let args = Args {
    ///     gxf: PathBuf::from("input.gtf"),
    ///     output: PathBuf::from("output.bed"),
    ///     // ... other fields
    ///     // Note: This is just a conceptual example
    /// };
    /// let config = Config::from_args(&args);
    /// ```
    pub fn from_args(args: &Args) -> Self {
        Self {
            input: args.gxf.clone(),
            output: args.output.clone(),
            threads: args.threads,
            parent_feature: args.parent_feature.clone(),
            child_features: args.child_features.clone(),
            parent_attribute: args.parent_attribute.clone(),
            child_attribute: args.child_attribute.clone(),
            bed_type: args.bed_type,
            additional_fields: args.additional_fields.clone(),
            chunks: args.chunks,
        }
    }
}
