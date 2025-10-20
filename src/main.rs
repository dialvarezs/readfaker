mod cli;

use clap::Parser;
use cli::Cli;
use readfaker::generator::{generate_reads, load_distributions};
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

    let (length_distribution, quality_distribution) = load_distributions(cli.input).unwrap();
    let reads = generate_reads(
        FastaReader::read(&cli.reference).unwrap(),
        length_distribution,
        quality_distribution,
        cli.num_reads,
        cli.seed,
    )
    .unwrap();

    let mut writer = FastqWriter::new(&cli.output).unwrap();
    writer.write_records(&reads).unwrap();
}
