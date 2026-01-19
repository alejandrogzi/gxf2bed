//! # gxf2bed
//!
//! The fastest GTF/GFF-to-BED converter chilling around.
//!
//! This library provides high-performance conversion of GTF and GFF3 files to various BED formats.
//! It leverages parallel processing and memory mapping for optimal performance on large genomic
//! annotation files.
//!
//! ## Usage
//!
//! ```rust, ignore
//! use gxf2bed::{Config, run};
//! use std::path::PathBuf;
//!
//! let config = Config {
//!     input: PathBuf::from("annotations.gtf"),
//!     output: PathBuf::from("output.bed"),
//!     threads: 4,
//!     parent_feature: Some("transcript".to_string()),
//!     child_features: Some(vec!["exon".to_string()]),
//!     parent_attribute: Some("transcript_id".to_string()),
//!     child_attribute: None,
//!     bed_type: gxf2bed::BedType::Bed12,
//!     additional_fields: None,
//!     chunks: 15000,
//! };
//!
//! let stats = run(&config)?;
//! println!("Conversion completed in {:?}", stats.elapsed);
//! println!("Memory used: {:.2} MB", stats.mem_delta_mb);
//! ```
//!
//! ## Examples
//!
//! ### Basic conversion
//!
//! ```rust, ignore
//! use gxf2bed::{Config, run, BedType};
//! use std::path::PathBuf;
//!
//! let config = Config {
//!     input: PathBuf::from("input.gtf"),
//!     output: PathBuf::from("output.bed"),
//!     threads: num_cpus::get(),
//!     parent_feature: Some("transcript".to_string()),
//!     child_features: Some(vec!["exon".to_string()]),
//!     parent_attribute: Some("transcript_id".to_string()),
//!     child_attribute: None,
//!     bed_type: BedType::Bed12,
//!     additional_fields: None,
//!     chunks: 15000,
//! };
//!
//! let stats = run(&config)?;
//! ```

//! ### Conversion with additional fields
//!
//! ```rust, ignore
//! use gxf2bed::{Config, run, BedType};
//! use std::path::PathBuf;
//!
//! let config = Config {
//!     input: PathBuf::from("input.gtf"),
//!     output: PathBuf::from("output.bed"),
//!     threads: 4,
//!     parent_feature: Some("transcript".to_string()),
//!     child_features: Some(vec!["exon".to_string()]),
//!     parent_attribute: Some("transcript_id".to_string()),
//!     child_attribute: None,
//!     bed_type: BedType::Bed12,
//!     additional_fields: Some(vec!["gene_name".to_string(), "gene_biotype".to_string()]),
//!     chunks: 15000,
//! };
//!
//! let stats = run(&config)?;
//! ```

//! ### Converting GFF3 files
//!
//! ```rust, ignore
//! use gxf2bed::{Config, run, BedType};
//! use std::path::PathBuf;
//!
//! let config = Config {
//!     input: PathBuf::from("input.gff3"),
//!     output: PathBuf::from("output.bed"),
//!     threads: 4,
//!     parent_feature: Some("mRNA".to_string()),
//!     child_features: Some(vec!["exon".to_string()]),
//!     parent_attribute: Some("ID".to_string()),
//!     child_attribute: None,
//!     bed_type: BedType::Bed12,
//!     additional_fields: None,
//!     chunks: 15000,
//! };
//!
//! let stats = run(&config)?;
//! ```

pub mod cli;
pub mod config;
pub mod convert;
pub mod detect;
pub mod error;
pub mod memory;

pub use cli::{Args, BedType};
pub use config::Config;
pub use convert::{run, RunStats};
pub use error::{Gxf2BedError, Result};
pub use memory::max_mem_usage_mb;
