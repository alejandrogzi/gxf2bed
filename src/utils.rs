use rayon::prelude::*;

use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use crate::gxf::GxfRecord;

pub fn reader<P: AsRef<Path> + Debug>(file: P) -> io::Result<String> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

pub fn parallel_parse<'a>(s: &'a str) -> Result<Vec<GxfRecord>, &'static str> {
    let records: Result<Vec<GxfRecord>, &'static str> =
        s.par_lines().map(|line| GxfRecord::parse(line)).collect();

    return records;
}

pub fn constructor<'a>(s: &'a str) -> Result<HashMap<String, HashMap<&str, String>>, &'static str> {
    s.par_lines()
        .map(|line| {
            if !line.starts_with("#") {
                Some(GxfRecord::parse(line))
            } else {
                None
            }
        })
        .filter_map(|x| x)
        .try_fold_with(HashMap::new(), |mut acc, record| {
            let record = record.unwrap();
            let tx_id = record.attr.transcript_id().to_owned();
            let entry = acc.entry(tx_id).or_insert(HashMap::new());

            if !record.chr.is_empty() {
                if record.feat == "transcript" {
                    entry.insert("chr", record.chr.to_owned());
                    entry.insert("start", record.start.to_string());
                    entry.insert("end", record.end.to_string());
                    entry.insert("strand", record.strand.to_string());
                } else if record.feat == "exon" {
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
        .unwrap_or(Err("Error"))
}

pub fn write_obj<P: AsRef<Path> + Debug>(filename: P, liner: Vec<(String, HashMap<&str, String>)>) {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
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
