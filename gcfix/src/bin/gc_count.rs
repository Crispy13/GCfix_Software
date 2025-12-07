use std::path::PathBuf;

use anyhow::Error;
use clap::Parser;
use gcfix::core::GCCounter;
use rayon::ThreadPoolBuilder;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// input bam path
    input_bam: String,

    /// out npy path
    out_npy: String,

    /// Minimum Mapping Quality for bam file read
    #[arg(short, long)]
    mapq: u8,

    /// minimum fragment length to include.
    #[arg(short, long)]
    start_len: i64,

    /// maximum fragment length to include.
    #[arg(short, long)]
    end_len: i64,

    /// Reference fasta file
    #[arg(short, long)]
    reference_fasta: String,

    /// Bin location csv file
    #[arg(short, long)]
    bin_location_csv: String,

    /// Number of threads.
    #[arg(short, long)]
    threads: usize,

    /// Lag
    #[arg(short, long, default_value_t = 10)]
    lag: usize,

    // #[command(subcommand)]
    // command: Option<Commands>,
}

fn main() -> Result<(), Error> {
    let cli = Cli::parse();

    let cgc = GCCounter::new(
        cli.mapq,
        cli.start_len,
        cli.end_len,
        cli.reference_fasta,
        cli.lag,
        cli.input_bam,
    );
    let tp = ThreadPoolBuilder::new().num_threads(cli.threads).build()?;

    
    tp.scope(|s| {

    });


    Ok(())
}
