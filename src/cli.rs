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

    /// Input FASTQ file to extract quality and length models
    #[arg(short = 'i', long, value_name = "FASTQ")]
    pub input: PathBuf,

    /// Output FASTQ file for simulated reads
    #[arg(short = 'o', long, value_name = "FASTQ")]
    pub output: PathBuf,

    /// Number of reads to generate
    #[arg(short = 'n', long, default_value = "100000")]
    pub num_reads: usize,

    /// Random seed for reproducibility
    #[arg(short = 's', long)]
    pub seed: Option<u64>,

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
        format!("{:<width$}", style(text).bold(), width = width + style(text).to_string().len() - text.len())
    }

    pub fn progress(text: impl Display) -> String {
        format!("{} {}", style("→").cyan(), style(text).dim())
    }

    pub fn success(text: impl Display) -> String {
        format!("{} {}", style("✓").green().bold(), style(text).green())
    }
}
