//! gxf2bed library module.

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
