use crate::cli::BedType;
use crate::config::Config;
use crate::detect::{detect_input_kind, Compression, InputFormat};
use crate::error::Result;
use crate::memory::max_mem_usage_mb;
use genepred::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed9, GenePred, Gff, Gtf, Reader, ReaderOptions};
use genepred::{Writer, WriterOptions};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{Duration, Instant};

/// Summary statistics for a conversion run.
#[derive(Debug, Clone, Copy)]
pub struct RunStats {
    /// Wall clock time spent in the conversion.
    pub elapsed: Duration,
    /// Delta in maximum RSS memory usage, in MB.
    pub mem_delta_mb: f64,
}

/// Runs a conversion with the provided configuration.
///
/// This function orchestrates the entire GTF/GFF to BED conversion process,
/// including input detection, parallel processing, and output writing.
///
/// # Arguments
///
/// * `config` - Configuration containing all conversion parameters
///
/// # Returns
///
/// Returns RunStats containing timing and memory usage information.
///
/// # Errors
///
/// Returns an error if any step of the conversion fails.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::{Config, run};
/// use std::path::PathBuf;
///
/// let config = Config {
///     input: PathBuf::from("input.gtf"),
///     output: PathBuf::from("output.bed"),
///     // ... other fields
/// };
/// let stats = run(&config)?;
/// println!("Conversion took: {:?}", stats.elapsed);
/// ```
pub fn run(config: &Config) -> Result<RunStats> {
    let start = Instant::now();
    let start_mem = max_mem_usage_mb();

    let reader_options = build_reader_options(config);
    let writer_options = build_writer_options(config);
    let input_kind = detect_input_kind(&config.input)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build()?;

    let chunks = pool.install(|| {
        process_input(
            &config.input,
            &input_kind,
            reader_options,
            config,
            &writer_options,
        )
    })?;

    write_output(&config.output, chunks)?;

    let elapsed = start.elapsed();
    let mem_delta = (max_mem_usage_mb() - start_mem).max(0.0);

    Ok(RunStats {
        elapsed,
        mem_delta_mb: mem_delta,
    })
}

/// Processes an input file based on the detected format.
///
/// This function dispatches to the appropriate reader type based on whether
/// the input is GTF or GFF format.
///
/// # Arguments
///
/// * `path` - Path to the input file
/// * `input_kind` - Detected format and compression information
/// * `reader_options` - Reader configuration options
/// * `config` - Main conversion configuration
/// * `writer_options` - Writer configuration options
///
/// # Returns
///
/// Returns processed data chunks ready for output.
///
/// # Errors
///
/// Returns an error if the input processing fails.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::process_input;
/// use std::path::Path;
/// // This function is typically called internally by run()
/// ```
fn process_input(
    path: &Path,
    input_kind: &crate::detect::InputKind,
    reader_options: ReaderOptions<'_>,
    config: &Config,
    writer_options: &WriterOptions,
) -> Result<Vec<(usize, Vec<u8>)>> {
    match input_kind.format {
        InputFormat::Gtf => process_reader::<Gtf>(
            path,
            reader_options,
            config,
            writer_options,
            input_kind.compression,
        ),
        InputFormat::Gff => process_reader::<Gff>(
            path,
            reader_options,
            config,
            writer_options,
            input_kind.compression,
        ),
    }
}

/// Processes a reader into in-memory BED chunks.
///
/// This function opens the appropriate reader, processes the input in parallel chunks,
/// and renders each chunk into BED format bytes.
///
/// # Arguments
///
/// * `path` - Path to the input file
/// * `reader_options` - Reader configuration options
/// * `config` - Main conversion configuration
/// * `writer_options` - Writer configuration options
/// * `compression` - Compression type of the input
///
/// # Returns
///
/// Returns ordered chunks of BED data ready for writing.
///
/// # Errors
///
/// Returns an error if reading or processing fails.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::process_reader;
/// use std::path::Path;
/// // This function is typically called internally by process_input()
/// ```
fn process_reader<F>(
    path: &Path,
    reader_options: ReaderOptions<'_>,
    config: &Config,
    writer_options: &WriterOptions,
    compression: Compression,
) -> Result<Vec<(usize, Vec<u8>)>>
where
    F: GxfReader + Send,
{
    let reader = open_reader::<F>(path, reader_options, compression)?;
    let mut outputs = reader
        .par_chunks(config.chunks)?
        .map(|(idx, chunk)| render_chunk(idx, chunk, config.bed_type, writer_options))
        .collect::<Vec<_>>();

    let mut merged = Vec::with_capacity(outputs.len());
    for output in outputs.drain(..) {
        merged.push(output?);
    }

    merged.sort_by_key(|(idx, _)| *idx);
    Ok(merged)
}

/// Opens a GTF/GFF reader with mmap or buffered I/O based on compression.
///
/// For compressed files, uses buffered I/O. For uncompressed files,
/// uses memory mapping for better performance.
///
/// # Arguments
///
/// * `path` - Path to the input file
/// * `options` - Reader configuration options
/// * `compression` - Compression type of the input
///
/// # Returns
///
/// Returns a configured reader for the specified format.
///
/// # Errors
///
/// Returns an error if the reader cannot be opened.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::open_reader;
/// use std::path::Path;
/// // This function is typically called internally by process_reader()
/// ```
fn open_reader<F: GxfReader>(
    path: &Path,
    options: ReaderOptions<'_>,
    compression: Compression,
) -> Result<Reader<F>> {
    if compression.is_compressed() {
        Ok(F::from_gxf_with_options(path, options)?)
    } else {
        Ok(F::from_mmap_with_options(path, options)?)
    }
}

/// Builds reader options from the config.
///
/// Constructs ReaderOptions based on the configuration parameters
/// for parent/child features and attributes.
///
/// # Arguments
///
/// * `config` - Configuration containing reader preferences
///
/// # Returns
///
/// Returns configured ReaderOptions for the genepred reader.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::build_reader_options;
/// use gxf2bed::Config;
/// let config = Config { /* ... */ };
/// let options = build_reader_options(&config);
/// ```
fn build_reader_options(config: &Config) -> ReaderOptions<'_> {
    let mut options = ReaderOptions::default();

    if let Some(parent) = &config.parent_feature {
        options = options.parent_feature(parent.as_bytes());
    }
    if let Some(childs) = &config.child_features {
        let filtered: Vec<&[u8]> = childs
            .iter()
            .filter(|c| !c.is_empty())
            .map(|c| c.as_bytes())
            .collect();
        if filtered.is_empty() {
            options = options.clear_child_features();
        } else {
            options = options.child_features(filtered);
        }
    }
    if let Some(parent) = &config.parent_attribute {
        options = options.parent_attribute(parent.as_bytes());
    }
    if let Some(child) = &config.child_attribute {
        options = options.child_attribute(child.as_bytes());
    }

    options
}

/// Builds writer options from the config.
///
/// Constructs WriterOptions based on additional field configuration.
/// If additional fields are specified, configures the allowlist
/// and enables non-numeric extras.
///
/// # Arguments
///
/// * `config` - Configuration containing writer preferences
///
/// # Returns
///
/// Returns configured WriterOptions for the BED writer.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::build_writer_options;
/// use gxf2bed::Config;
/// let config = Config { /* ... */ };
/// let options = build_writer_options(&config);
/// ```
fn build_writer_options(config: &Config) -> WriterOptions {
    if let Some(additional_fields) = &config.additional_fields {
        let mut options = WriterOptions::default();
        options = options.extras_allowlist(additional_fields.clone());
        options = options.include_non_numeric_extras(true);
        options
    } else {
        WriterOptions::default()
    }
}

/// Renders a chunk of parsed records into BED bytes.
///
/// Converts a chunk of GenePred records into the specified BED format
/// and returns the rendered bytes along with the chunk index.
///
/// # Arguments
///
/// * `idx` - Index of this chunk for ordering
/// * `chunk` - Vector of parsed GenePred records
/// * `bed_type` - Target BED format type
/// * `writer_options` - Writer configuration options
///
/// # Returns
///
/// Returns a tuple containing the chunk index and rendered BED bytes.
///
/// # Errors
///
/// Returns an error if record parsing or writing fails.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::render_chunk;
/// // This function is typically called internally by process_reader()
/// ```
fn render_chunk(
    idx: usize,
    chunk: Vec<genepred::ReaderResult<GenePred>>,
    bed_type: BedType,
    writer_options: &WriterOptions,
) -> Result<(usize, Vec<u8>)> {
    let mut buffer = Vec::with_capacity(chunk.len().saturating_mul(128));
    {
        let mut writer = BufWriter::with_capacity(128 * 1024, &mut buffer);
        for record in chunk {
            let record = record?;
            write_record(&record, &mut writer, bed_type, writer_options)?;
        }
        writer.flush()?;
    }
    Ok((idx, buffer))
}

/// Writes a single record in the configured BED format.
///
/// Converts a GenePred record to the specified BED format and writes it
/// to the provided writer using the configured options.
///
/// # Arguments
///
/// * `record` - GenePred record to convert
/// * `writer` - Output writer to write the BED record to
/// * `bed_type` - Target BED format type
/// * `writer_options` - Writer configuration options
///
/// # Returns
///
/// Returns Ok(()) on successful write.
///
/// # Errors
///
/// Returns an error if the record cannot be written in the specified format.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::write_record;
/// use std::io::Write;
/// // This function is typically called internally by render_chunk()
/// ```
fn write_record<W: Write>(
    record: &GenePred,
    writer: &mut W,
    bed_type: BedType,
    writer_options: &WriterOptions,
) -> Result<()> {
    match bed_type {
        BedType::Bed3 => Writer::<Bed3>::from_record_with_options(record, writer, writer_options)?,
        BedType::Bed4 => Writer::<Bed4>::from_record_with_options(record, writer, writer_options)?,
        BedType::Bed5 => Writer::<Bed5>::from_record_with_options(record, writer, writer_options)?,
        BedType::Bed6 => Writer::<Bed6>::from_record_with_options(record, writer, writer_options)?,
        BedType::Bed9 => Writer::<Bed9>::from_record_with_options(record, writer, writer_options)?,
        BedType::Bed12 => {
            Writer::<Bed12>::from_record_with_options(record, writer, writer_options)?
        }
    }
    Ok(())
}

/// Writes ordered chunks to the output path.
///
/// Creates the output file and writes all chunks in order,
/// combining the parallel processed data into a single BED file.
///
/// # Arguments
///
/// * `path` - Output file path to write to
/// * `chunks` - Ordered vector of (index, data) tuples
///
/// # Returns
///
/// Returns Ok(()) on successful write.
///
/// # Errors
///
/// Returns an error if file creation or writing fails.
///
/// # Example
///
/// ```rust, ignore
/// use gxf2bed::convert::write_output;
/// use std::path::Path;
/// let chunks = vec![(0, b"chr1\t100\t200\n".to_vec())];
/// write_output(Path::new("output.bed"), chunks)?;
/// ```
fn write_output(path: &Path, chunks: Vec<(usize, Vec<u8>)>) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = BufWriter::with_capacity(256 * 1024, file);
    for (_, buffer) in chunks {
        writer.write_all(&buffer)?;
    }
    writer.flush()?;
    Ok(())
}

/// Trait for opening GTF/GFF readers with custom options.
///
/// This trait abstracts over GTF and GFF readers, providing a common interface
/// for opening readers with custom configuration options.
trait GxfReader: genepred::BedFormat + Into<GenePred> + Sized {
    /// Opens a buffered reader with custom options.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the input file
    /// * `options` - Reader configuration options
    ///
    /// # Returns
    ///
    /// Returns a Reader for the specified format.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::convert::GxfReader;
    /// use genepred::Gtf;
    /// let reader = Gtf::from_gxf_with_options("file.gtf", options)?;
    /// ```
    fn from_gxf_with_options<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> genepred::ReaderResult<Reader<Self>>;

    /// Opens a memory-mapped reader with custom options.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the input file
    /// * `options` - Reader configuration options
    ///
    /// # Returns
    ///
    /// Returns a memory-mapped Reader for the specified format.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::convert::GxfReader;
    /// use genepred::Gtf;
    /// let reader = Gtf::from_mmap_with_options("file.gtf", options)?;
    /// ```
    fn from_mmap_with_options<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> genepred::ReaderResult<Reader<Self>>;
}

impl GxfReader for Gtf {
    /// Opens a buffered GTF reader with custom options.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the GTF file
    /// * `options` - Reader configuration options
    ///
    /// # Returns
    ///
    /// Returns a buffered GTF Reader.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::convert::GxfReader;
    /// use genepred::Gtf;
    /// let reader = Gtf::from_gxf_with_options("file.gtf", options)?;
    /// ```
    fn from_gxf_with_options<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> genepred::ReaderResult<Reader<Self>> {
        Reader::<Gtf>::from_gxf_with_options(path, options)
    }

    /// Opens a memory-mapped GTF reader with custom options.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the GTF file
    /// * `options` - Reader configuration options
    ///
    /// # Returns
    ///
    /// Returns a memory-mapped GTF Reader.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::convert::GxfReader;
    /// use genepred::Gtf;
    /// let reader = Gtf::from_mmap_with_options("file.gtf", options)?;
    /// ```
    fn from_mmap_with_options<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> genepred::ReaderResult<Reader<Self>> {
        Reader::<Gtf>::from_mmap_with_options(path, options)
    }
}

impl GxfReader for Gff {
    /// Opens a buffered GFF reader with custom options.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the GFF file
    /// * `options` - Reader configuration options
    ///
    /// # Returns
    ///
    /// Returns a buffered GFF Reader.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::convert::GxfReader;
    /// use genepred::Gff;
    /// let reader = Gff::from_gxf_with_options("file.gff", options)?;
    /// ```
    fn from_gxf_with_options<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> genepred::ReaderResult<Reader<Self>> {
        Reader::<Gff>::from_gxf_with_options(path, options)
    }

    /// Opens a memory-mapped GFF reader with custom options.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the GFF file
    /// * `options` - Reader configuration options
    ///
    /// # Returns
    ///
    /// Returns a memory-mapped GFF Reader.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use gxf2bed::convert::GxfReader;
    /// use genepred::Gff;
    /// let reader = Gff::from_mmap_with_options("file.gff", options)?;
    /// ```
    fn from_mmap_with_options<P: AsRef<Path>>(
        path: P,
        options: ReaderOptions<'_>,
    ) -> genepred::ReaderResult<Reader<Self>> {
        Reader::<Gff>::from_mmap_with_options(path, options)
    }
}
