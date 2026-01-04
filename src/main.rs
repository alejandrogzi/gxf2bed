//! # gxf2bed
//!
//! The fastest GTF/GFF-to-BED converter chilling around.
//!
//! ## Features
//!
//! - Converts GTF and GFF3 files to BED format
//! - Supports multiple BED formats (BED3, BED4, BED5, BED6, BED9, BED12)
//! - Handles compressed input files (gzip, zstd, bzip2)
//! - Multi-threaded processing for maximum performance
//! - Customizable parent/child feature and attribute mapping
//! - Optional additional fields from GTF/GFF attributes
//! - Memory-mapped I/O for uncompressed files
//!
//! ## Usage
//!
//! ```bash
//! gxf2bed -i <INPUT> -o <OUTPUT> [OPTIONS]
//!
//! Required arguments:
//!   -i, --input <GXF>          Path to GTF/GFF file
//!   -o, --output <BED>         Path to output BED file
//!
//! Optional arguments:
//!   -T, --threads <THREADS>    Number of threads (default: CPU count)
//!   -F, --parent-feature <PARENT>     Parent feature
//!   -f, --child-features <CHILDS>    Child features (comma-separated)
//!   -A, --parent-attribute <FEATURE>  Feature to extract
//!   -a, --child-attribute <CHILD>     Child feature to extract
//!   -t, --type <BED_TYPE>      BED type format (3, 4, 5, 6, 9, 12) [default: 12]
//!   -d, --additional-fields <ADDITIONAL>  BED additional fields (comma-separated)
//!   -c, --chunks <CHUNKS>      Chunk size for parallel processing [default: 15000]
//!   -h, --help                 Print help
//!   -V, --version              Print version
//! ```
//!
//! ## Examples
//!
//! ### Basic GTF to BED conversion
//!
//! ```bash
//! gxf2bed -i annotations.gtf -o output.bed
//! ```
//!
//! ### Convert with custom threads and BED format
//!
//! ```bash
//! gxf2bed -i annotations.gtf -o output.bed -T 8 -t 6
//! ```
//!
//! ### Convert compressed GFF3 file
//!
//! ```bash
//! gxf2bed -i annotations.gff3.gz -o output.bed
//! ```
//!
//! ### Custom parent/child features for GFF3
//!
//! ```bash
//! gxf2bed -i annotations.gff3 -o output.bed -F mRNA -f exon -A ID
//! ```
//!
//! ### Include additional fields from GTF attributes
//!
//! ```bash
//! gxf2bed -i annotations.gtf -o output.bed -d gene_name,gene_biotype
//! ```
//!
//! ```bash
//! gxf2bed -i annotations.gff -o output.bed -F "" -f ""
//! ```
use gxf2bed::{run, Args, Config};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    simple_logger::init_with_level(Level::Info).unwrap();

    let args = Args::parse();
    log::info!("{:?}", args);

    let config = Config::from_args(&args);
    log::info!("Using {} threads", config.threads);

    let stats = run(&config)?;
    log::info!("Elapsed: {:.4?} secs", stats.elapsed.as_secs_f32());
    log::info!("Memory: {:.2} MB", stats.mem_delta_mb);

    Ok(())
}
