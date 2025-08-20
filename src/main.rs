//! The fastest GTF/GFF-to-BED converter chilling around
//! Alejandro Gonzales-Irribarren, 2025

use clap::Parser;
use log::Level;

use gxf2bed::{
    cli::Args,
    utils::{convert, initialize},
};

fn main() {
    initialize();
    let st = std::time::Instant::now();
    simple_logger::init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();
    args.check().unwrap_or_else(|e| {
        log::error!("{}", e);
        std::process::exit(1);
    });

    log::info!("{:?}", args);

    convert(args);

    log::info!("Elapsed: {:.4?} secs", st.elapsed().as_secs_f32());
}
