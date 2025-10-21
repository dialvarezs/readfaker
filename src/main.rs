mod cli;

use clap::Parser;
use cli::Cli;
use readfaker::generator::{load_distributions, ReadGenerator};
use readfaker::io::{FastaReader, FastqWriter};

fn main() {
    let cli = Cli::parse();

    if cli.verbose {
        eprintln!("=== ReadFaker Configuration ===");
        eprintln!("Reference: {}", cli.reference.display());
        eprintln!("Input: {}", cli.input.display());
        eprintln!("Output: {}", cli.output.display());
        eprintln!("Number of reads: {}", cli.num_reads);
        if let Some(seed) = cli.seed {
            eprintln!("Random seed: {}", seed);
        }
        eprintln!("================================\n");
    }

    let (length_distribution, quality_distribution) = load_distributions(&cli.input).unwrap();
    let mut generator = ReadGenerator::new(
        FastaReader::read(&cli.reference).unwrap(),
        length_distribution,
        quality_distribution,
        cli.seed,
    );
    let mut writer = FastqWriter::new(&cli.output).unwrap();

    for _ in 0..cli.num_reads {
        let read = generator.generate_read().unwrap();
        writer.write_record(&read).unwrap();
    }
}
