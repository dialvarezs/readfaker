use anyhow::Result;
use clap::Parser;
use readfaker::cli::{fmt, Cli};
use readfaker::generator::ReadGenerator;
use readfaker::io::{FastaReader, FastqWriter};
use readfaker::models::ErrorModel;
use readfaker::utils::load_models;

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.verbose {
        eprintln!("{}", fmt::header("ReadFaker Configuration"));
        eprintln!("{}: {}", fmt::param_aligned("Reference", 16), cli.reference.display());
        eprintln!("{}: {}", fmt::param_aligned("Input", 16), cli.input.display());
        eprintln!("{}: {}", fmt::param_aligned("Output", 16), cli.output.display());
        eprintln!("{}: {}", fmt::param_aligned("Number of reads", 16), cli.num_reads);
        if let Some(seed) = cli.seed {
            eprintln!("{}: {}", fmt::param_aligned("Random seed", 16), seed);
        }
        eprintln!();
    }

    if cli.verbose {
        eprintln!("{}", fmt::progress("Creating models from input FASTQ..."));
    }
    let (length_model, quality_model) = load_models(&cli.input, cli.seed)?;

    let error_model = ErrorModel::new(None, None, None)?;
    let mut generator = ReadGenerator::new(
        FastaReader::read(&cli.reference)?,
        length_model,
        quality_model,
        error_model,
        cli.seed,
    )?;
    let mut writer = FastqWriter::new(&cli.output)?;

    if cli.verbose {
        eprintln!("{}", fmt::progress(format!("Generating {} reads...", cli.num_reads)));
    }

    for _ in 0..cli.num_reads {
        let read = generator.generate_read()?;
        writer.write_record(&read)?;
    }

    if cli.verbose {
        eprintln!("{}", fmt::success(format!("Output written to {}", cli.output.display())));
    }

    Ok(())
}
