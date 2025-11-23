use anyhow::Result;
use clap::Parser;
use readfaker::cli::{Cli, fmt};
use readfaker::generator::ReadGenerator;
use readfaker::io::{BamWriter, FastaReader, FastqWriter};
use readfaker::models::ErrorModel;
use readfaker::utils::load_models;

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.verbose {
        eprintln!("{}", fmt::header("ReadFaker Configuration"));
        eprintln!(
            "{}: {}",
            fmt::param_aligned("Reference", 16),
            cli.reference.display()
        );
        eprintln!(
            "{}: {}",
            fmt::param_aligned("Input", 16),
            cli.input.display()
        );
        eprintln!(
            "{}: {}",
            fmt::param_aligned("Output", 16),
            cli.output.display()
        );
        eprintln!(
            "{}: {}",
            fmt::param_aligned("Number of reads", 16),
            cli.num_reads
        );
        if let Some(seed) = cli.seed {
            eprintln!("{}: {}", fmt::param_aligned("Random seed", 16), seed);
        }
        eprintln!();
    }

    if cli.verbose {
        eprintln!("{}", fmt::progress("Creating models from input FASTQ..."));
    }
    let (length_model, quality_model) = load_models(&cli.input, cli.seed)?;

    let error_model = ErrorModel::new(
        cli.error_sub,
        cli.error_ins,
        cli.error_del,
        cli.error_ins_ext,
        cli.error_del_ext,
    )?;

    if cli.verbose {
        eprintln!("Error Model Configuration:");
        eprintln!(
            "{}: {:.2}",
            fmt::param_aligned("Substitution rate", 20),
            error_model.substitution_rate
        );
        eprintln!(
            "{}: {:.2}",
            fmt::param_aligned("Insertion rate", 20),
            error_model.insertion_rate
        );
        eprintln!(
            "{}: {:.2}",
            fmt::param_aligned("Deletion rate", 20),
            error_model.deletion_rate
        );
        eprintln!(
            "{}: {:.2}",
            fmt::param_aligned("Ins. extension rate", 20),
            error_model.insertion_extension_rate
        );
        eprintln!(
            "{}: {:.2}",
            fmt::param_aligned("Del. extension rate", 20),
            error_model.deletion_extension_rate
        );
        eprintln!();
    }

    let mut generator = ReadGenerator::new(
        FastaReader::read(&cli.reference)?,
        length_model,
        quality_model,
        error_model,
        cli.seed,
    )?;

    // Detect output format based on extension
    let output_ext = cli
        .output
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("");

    if cli.verbose {
        eprintln!(
            "{}",
            fmt::progress(format!("Generating {} reads...", cli.num_reads))
        );
    }

    match output_ext.to_lowercase().as_str() {
        "bam" => {
            let mut writer = BamWriter::new(&cli.output)?;
            for _ in 0..cli.num_reads {
                let read = generator.generate_read()?;
                let name = std::str::from_utf8(read.name()).unwrap_or("unknown");
                writer.write_record(name, read.sequence(), read.quality_scores())?;
            }
            writer.finish()?;
        }
        _ => {
            // Default to FASTQ for all other extensions
            let mut writer = FastqWriter::new(&cli.output, cli.compression_threads)?;
            for _ in 0..cli.num_reads {
                let read = generator.generate_read()?;
                writer.write_record(&read)?;
            }
            writer.finish()?;
        }
    }

    if cli.verbose {
        eprintln!(
            "{}",
            fmt::success(format!("Output written to {}", cli.output.display()))
        );
    }

    Ok(())
}
