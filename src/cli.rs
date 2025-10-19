// Command-line interface definition

use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "readfaker",
    version,
    about = "Simulate Oxford Nanopore reads with realistic quality profiles",
    long_about = None
)]
pub struct Cli {
    /// Reference sequences (FASTA format) to sample reads from
    #[arg(short = 'r', long, value_name = "FASTA")]
    pub reference: PathBuf,

    /// Input FASTQ file to extract quality and length distributions
    #[arg(short = 'i', long, value_name = "FASTQ")]
    pub input: PathBuf,

    /// Output FASTQ file for simulated reads
    #[arg(short = 'o', long, value_name = "FASTQ")]
    pub output: PathBuf,

    /// Number of reads to generate
    #[arg(short = 'n', long, default_value = "10000")]
    pub num_reads: usize,

    /// Random seed for reproducibility
    #[arg(short = 's', long)]
    pub seed: Option<u64>,

    /// Enable verbose output
    #[arg(short, long)]
    pub verbose: bool,
}
