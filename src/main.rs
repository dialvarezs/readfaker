use anyhow::Result;
use clap::Parser;
use readfaker::cli::{Cli, fmt};
use readfaker::generator::{ORIGIN_TAG, ReadGenerator};
use readfaker::io::{BamWriter, CompressionLevel, FastaReader, FastqWriter};
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
        for (name, rate) in [
            ("Substitution rate", error_model.substitution_rate),
            ("Insertion rate", error_model.insertion_rate),
            ("Deletion rate", error_model.deletion_rate),
            ("Ins. extension rate", error_model.insertion_extension_rate),
            ("Del. extension rate", error_model.deletion_extension_rate),
        ] {
            eprintln!("{}: {:.2}", fmt::param_aligned(name, 20), rate);
        }
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
        .unwrap_or("")
        .to_lowercase();

    if cli.verbose {
        eprintln!(
            "{}",
            fmt::progress(format!("Generating {} reads...", cli.num_reads))
        );
    }

    let compression_level = CompressionLevel::try_from(cli.compression_level)?;
    match output_ext.as_str() {
        "bam" => {
            let origin_prefix = format!("{ORIGIN_TAG}:Z:");
            let mut writer =
                BamWriter::new(&cli.output, cli.compression_threads, compression_level)?;
            for _ in 0..cli.num_reads {
                let read = generator.generate_read()?;
                let name = std::str::from_utf8(read.name()).expect("UUID should be valid UTF-8");
                let origin = read
                    .description()
                    .strip_prefix(origin_prefix.as_bytes())
                    .filter(|_| !cli.no_origin)
                    .map(|o| std::str::from_utf8(o).expect("origin should be valid UTF-8"));
                writer.write_record(name, read.sequence(), read.quality_scores(), origin)?;
            }
            writer.finish()?;
        }
        _ => {
            // Default to FASTQ for all other extensions
            let mut writer =
                FastqWriter::new(&cli.output, cli.compression_threads, compression_level)?;
            for _ in 0..cli.num_reads {
                let mut read = generator.generate_read()?;
                if cli.no_origin {
                    read.description_mut().clear();
                }
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
