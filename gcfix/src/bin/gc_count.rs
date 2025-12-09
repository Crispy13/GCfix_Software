use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    ops::AddAssign,
    path::PathBuf,
    str::{FromStr, Split},
};

use anyhow::{Error, anyhow};
use clap::Parser;
use crackle_kit::tracing::{self, Level, event};
use crackle_kit::{tracing::level_filters::LevelFilter, tracing_kit::setup_logging_stderr_only};
use gcfix::core::GCCounter;
use ndarray::Array2;
use numpy::IntoPyArray;
use pyo3::{Python, types::PyAnyMethods};
use rayon::{
    ThreadPoolBuilder,
    iter::{
        IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
    },
};

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
    
    /// Log level, Default: INFO
    #[arg(long, default_value_t = LevelFilter::INFO)]
    log_level: LevelFilter,

    // #[command(subcommand)]
    // command: Option<Commands>,
}

fn main() -> Result<(), Error> {
    let cli = Cli::parse();

    setup_logging_stderr_only(cli.log_level)?;


    let cgc = GCCounter::new(
        cli.mapq,
        cli.start_len,
        cli.end_len,
        cli.reference_fasta,
        cli.lag,
        cli.input_bam,
    )?;

    // build bin location vector
    let bin_location_vec = {
        let mut res = Vec::with_capacity(48000);
        let mut br = BufReader::new(File::open(&cli.bin_location_csv)?);
        let mut buf = String::new();

        // read first line. it's header.
        br.read_line(&mut buf)?;

        macro_rules! parse_pos {
            ($spl:ident) => {{
                let s = $spl.next().unwrap();
                match s.parse::<i64>() {
                    Ok(v) => Ok(v),
                    Err(err) => Err(Error::from(err).context(s.to_owned())),
                }
            }};
        }

        while let Ok(true) = {
            buf.clear();
            br.read_line(&mut buf).map(|n| n > 0)
        } {
            let mut spl = buf.trim().split(",");
            res.push((
                spl.next().unwrap().to_owned(),
                parse_pos!(spl)?,
                parse_pos!(spl)?,
            ));
        }

        res
    };

    let tp = ThreadPoolBuilder::new().num_threads(cli.threads).build()?;

    let arr_iden_fn = || {
        Array2::zeros(GCCounter::res_arr_shape(
            cgc.start_len as usize,
            cgc.end_len as usize,
        ))
    };

    let chunk_size = (bin_location_vec.len() / cli.threads).min(64).max(1);

    debug_assert!(chunk_size > 0);

    event!(Level::DEBUG, "chunk_size={chunk_size}");
    let res_arr = tp.scope(|_s| {
        let r = bin_location_vec
            .par_iter()
            .chunks(chunk_size)
            .map(|b| {
                let mut res_arr = arr_iden_fn();
                b.into_iter().try_for_each(|r| {
                    res_arr += &cgc.count_gc(r)?;
                    Ok::<_, Error>(())
                })?;

                Ok::<_, Error>(res_arr)
            })
            .try_reduce(arr_iden_fn, |mut a, b| {
                a.add_assign(&b);
                Ok(a)
            })?;

        Ok::<_, Error>(r)
    })?;

    Python::attach(|py| {
        let np = py.import("numpy")?;
        let py_arr = res_arr.into_pyarray(py);

        np.getattr("save")?.call1((cli.out_npy, py_arr))?;

        Ok::<_, Error>(())
    })?;

    Ok(())
}
