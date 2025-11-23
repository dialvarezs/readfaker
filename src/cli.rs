use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use console::style;
use std::fmt::Display;
use std::path::PathBuf;

fn get_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Cyan.on_default() | Effects::BOLD)
        .usage(AnsiColor::Cyan.on_default() | Effects::BOLD)
        .literal(AnsiColor::Green.on_default() | Effects::BOLD)
        .placeholder(AnsiColor::Yellow.on_default())
}

#[derive(Parser)]
#[command(
    name = "readfaker",
    version,
    about = "Simulate Oxford Nanopore reads with realistic quality profiles",
    long_about = None,
    styles = get_styles()
)]
pub struct Cli {
    /// Reference sequences (FASTA format) to sample reads from
    #[arg(short = 'r', long, value_name = "FASTA")]
    pub reference: PathBuf,

    /// Input file to extract quality and length models (FASTQ or BAM)
    #[arg(short = 'i', long, value_name = "FILE")]
    pub input: PathBuf,

    /// Output file for simulated reads (FASTQ or BAM, detected by extension)
    #[arg(short = 'o', long, value_name = "FILE")]
    pub output: PathBuf,

    /// Number of reads to generate
    #[arg(short = 'n', long, default_value = "100000")]
    pub num_reads: usize,

    /// Random seed for reproducibility
    #[arg(short = 's', long)]
    pub seed: Option<u64>,

    /// Number of compression threads (default: 4)
    #[arg(long = "compression-threads", default_value = "4")]
    pub compression_threads: Option<usize>,

    /// Error substitution rate (default: 0.7)
    #[arg(long, value_name = "RATE")]
    pub error_sub: Option<f64>,

    /// Error insertion rate (default: 0.1)
    #[arg(long, value_name = "RATE")]
    pub error_ins: Option<f64>,

    /// Error deletion rate (default: 0.2)
    #[arg(long, value_name = "RATE")]
    pub error_del: Option<f64>,

    /// Error insertion extension rate (default: 0.4)
    #[arg(long, value_name = "RATE")]
    pub error_ins_ext: Option<f64>,

    /// Error deletion extension rate (default: 0.4)
    #[arg(long, value_name = "RATE")]
    pub error_del_ext: Option<f64>,

    /// Enable verbose output
    #[arg(short, long)]
    pub verbose: bool,
}

/// Formatting utilities for console output
pub mod fmt {
    use super::*;

    pub fn header(text: impl Display) -> String {
        style(text).bold().cyan().to_string()
    }

    pub fn param(text: impl Display) -> String {
        style(text).bold().to_string()
    }

    /// Format a parameter with proper alignment (width must account for raw text length)
    pub fn param_aligned(text: &str, width: usize) -> String {
        format!(
            "{:<width$}",
            style(text).bold(),
            width = width + style(text).to_string().len() - text.len()
        )
    }

    pub fn progress(text: impl Display) -> String {
        format!("{} {}", style("→").cyan(), style(text).dim())
    }

    pub fn success(text: impl Display) -> String {
        format!("{} {}", style("✓").green().bold(), style(text).green())
    }
}
