mod cli;

use clap::Parser;
use cli::Cli;
use readfaker::generator::load_distributions;

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

    load_distributions(cli.input).unwrap();
}
