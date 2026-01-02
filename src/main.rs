use clap::Parser;
use log::Level;

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
