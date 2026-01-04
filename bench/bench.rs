use clap::Parser;
use std::path::PathBuf;
use std::process::{Command, ExitStatus, Stdio};

const TOOLS: [&str; 2] = [
    "target/release/gxf2bed -i {bed} -o ./test/output.bed",
    "target/release/gxf2bed -i {bed} -o ./test/output.bed.gz",
];
const GTFILES: [&str; 4] = [
    "Homo_sapiens.GRCh38.112.gtf",
    "Homo_sapiens.GRCh38.112.gtf.gz",
    "gencode.v46.chr_patch_hapl_scaff.annotation.gff3.gz",
    "gencode.v46.chr_patch_hapl_scaff.annotation.gff3",
];
const STDOUT: &str = "*.bed";
const CSV: &str = "bed.files.csv";
const MD: &str = "bed.files.md";

#[derive(Debug, Parser)]
pub struct Args {
    #[clap(
        short = 'd',
        long = "dir",
        help = "Path to the test files directory",
        default_value = "assets"
    )]
    assets: PathBuf,

    #[clap(short = 'a',
        value_delimiter = ',',
        num_args = 1..,
        help = "Extra arguments to pass to hyperfine"
    )]
    hyperfine_args: Vec<String>,
}

pub struct HyperfineCall {
    pub warmup: u32,
    pub min_runs: u32,
    pub max_runs: Option<u32>,
    pub export_csv: Option<String>,
    pub export_markdown: Option<String>,
    pub parameters: Vec<(String, Vec<String>)>,
    pub setup: Option<String>,
    pub cleanup: Option<String>,
    pub commands: Vec<String>,
    pub extras: Vec<String>,
}

impl Default for HyperfineCall {
    fn default() -> Self {
        Self {
            warmup: 3,
            min_runs: 5,
            max_runs: None,
            export_csv: None,
            export_markdown: None,
            parameters: Vec::new(),
            setup: None,
            cleanup: None,
            commands: Vec::new(),
            extras: Vec::new(),
        }
    }
}

impl HyperfineCall {
    pub fn invoke(&self) -> ExitStatus {
        let mut command = Command::new("hyperfine");

        command
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit())
            .stdin(Stdio::null());

        command.arg("--warmup").arg(self.warmup.to_string());
        command.arg("--min-runs").arg(self.min_runs.to_string());
        if let Some(export_csv) = &self.export_csv {
            command.arg("--export-csv").arg(export_csv);
        }
        if let Some(export_markdown) = &self.export_markdown {
            command.arg("--export-markdown").arg(export_markdown);
        }
        for (flag, values) in &self.parameters {
            command.arg("-L").arg(flag).arg(values.join(","));
        }
        if let Some(setup) = &self.setup {
            command.arg("--setup").arg(setup);
        }
        if let Some(cleanup) = &self.cleanup {
            command.arg("--cleanup").arg(cleanup);
        }
        if let Some(max_runs) = self.max_runs {
            command.arg("--max-runs").arg(max_runs.to_string());
        }
        if !self.extras.is_empty() {
            command.args(&self.extras);
        }

        for cmd in &self.commands {
            command.arg(cmd);
        }

        command.status().expect("Failed to run hyperfine")
    }
}

fn benchmark() -> Result<(String, String), Box<dyn std::error::Error>> {
    let args = Args::parse();

    std::fs::create_dir_all("runs")?;
    let assets = args.assets.to_string_lossy();

    #[allow(clippy::needless_update)]
    let code = HyperfineCall {
        warmup: 5,
        min_runs: 5,
        max_runs: Some(20),
        export_csv: Some(format!("runs/{}", CSV).to_string()),
        export_markdown: Some(format!("runs/{}", MD).to_string()),
        parameters: vec![(
            "bed".to_string(),
            GTFILES
                .iter()
                .map(|s| format!("{}/{}", assets, s))
                .collect(),
        )],
        setup: Some("cargo build --release".to_string()),
        cleanup: Some(format!("rm -f {} test/*.bed", STDOUT)),
        commands: TOOLS
            .iter()
            .map(|cmd| cmd.to_string())
            .collect::<Vec<String>>(),
        extras: args
            .hyperfine_args
            .iter()
            .map(|s| format!("--{}", s))
            .collect(),
        ..Default::default()
    }
    .invoke()
    .code()
    .expect("Benchmark terminated unexpectedly");

    if code != 0 {
        return Err(format!("Benchmark failed with exit code {}", code).into());
    }

    Ok((format!("runs/{}", CSV), format!("runs/{}", MD)))
}

fn main() {
    match benchmark() {
        Ok((csv, md)) => {
            println!("Benchmark results saved to:");
            println!("  - {}", csv);
            println!("  - {}", md);
        }
        Err(e) => eprintln!("{}", e),
    }
}
