use std::error::Error;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use hashbrown::HashMap;

use flate2::{read::GzDecoder, write::GzEncoder, Compression};

use colored::Colorize;
use indoc::indoc;

const VERSION: &str = env!("CARGO_PKG_VERSION");

pub fn raw<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let mut file = File::open(f)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    Ok(contents)
}

pub fn with_gz<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let file = File::open(f)?;
    let mut decoder = GzDecoder::new(file);
    let mut contents = String::new();

    decoder.read_to_string(&mut contents)?;
    Ok(contents)
}

pub fn write_obj<P: AsRef<Path> + Debug>(filename: P, liner: Vec<(String, HashMap<&str, String>)>) {
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

    let mut warn_missing_child_count = 0;
    for x in liner.iter() {
        let start = x.1.get("start").unwrap().parse::<u32>().unwrap();
        let start_codon = x.1.get("start_codon").unwrap_or(x.1.get("start").unwrap());
        let end_codon = x.1.get("stop_codon").unwrap_or(x.1.get("end").unwrap());

        let mut exon_sizes = if let Some(exon_sizes) = x.1.get("exon_sizes") {
            exon_sizes
                .split(',')
                .filter(|x| !x.is_empty())
                .collect::<Vec<&str>>()
        } else {
            warn_missing_child_count += 1;
            Vec::new()
        };

        let mut exon_starts = if let Some(exon_starts) = x.1.get("exon_starts") {
            exon_starts
                .split(',')
                .filter(|x| !x.is_empty())
                .map(|x| x.parse::<u32>().unwrap() - start)
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
        } else {
            warn_missing_child_count += 1;
            Vec::new()
        };

        if exon_sizes.is_empty() || exon_starts.is_empty() {
            continue;
        }

        if x.1.get("strand") == Some(&String::from("-")) {
            exon_sizes.reverse();
            exon_starts.reverse();
        }

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            x.1["chr"],
            x.1["start"],
            x.1["end"],
            x.0,
            "0",
            x.1["strand"],
            start_codon,
            end_codon,
            "0",
            x.1["exons"].len(),
            exon_sizes.join(",") + ",",
            exon_starts.join(",") + ",",
        );
        writeln!(writer, "{}", line).unwrap();
    }

    log::warn!(
        "{} entries were skipped due to missing child!",
        warn_missing_child_count
    );

    writer.flush().unwrap();
    log::info!("Done writing!");
}

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

pub fn msg() {
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
