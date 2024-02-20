// use rayon::prelude::*;

use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use colored::Colorize;
use indoc::indoc;

// use crate::gxf::GxfRecord;

const VERSION: &str = env!("CARGO_PKG_VERSION");

pub fn reader<P: AsRef<Path> + Debug>(file: P) -> io::Result<String> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

// pub fn parallel_parse<'a>(s: &'a str) -> Result<Vec<GxfRecord>, &'static str> {
//     let records: Result<Vec<GxfRecord>, &'static str> =
//         s.par_lines().map(|line| GxfRecord::parse(line)).collect();
//
//     return records;
// }

pub fn write_obj<P: AsRef<Path> + Debug>(filename: P, liner: Vec<(String, HashMap<&str, String>)>) {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    log::info!("Writing to {:?}", filename);
    let mut writer = BufWriter::new(f);

    for x in liner.iter() {
        let start = x.1.get("start").unwrap().parse::<u32>().unwrap();
        let start_codon = x.1.get("start_codon").unwrap_or(x.1.get("start").unwrap());
        let end_codon = x.1.get("stop_codon").unwrap_or(x.1.get("end").unwrap());

        let mut exon_sizes =
            x.1.get("exon_sizes")
                .unwrap()
                .split(",")
                .filter(|x| !x.is_empty())
                .collect::<Vec<&str>>();

        let mut exon_starts =
            x.1.get("exon_starts")
                .unwrap()
                .split(",")
                .filter(|x| !x.is_empty())
                .map(|x| x.parse::<u32>().unwrap() - start)
                .map(|x| x.to_string())
                .collect::<Vec<String>>();

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
