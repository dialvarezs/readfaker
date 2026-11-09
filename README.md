# ReadFaker

A simple and fast Oxford Nanopore sequencing reads with realistic quality profiles by extracting empirical models from real sequencing data.

## Features

- Creates empirical models for read length and quality scores.
- Supports both FASTQ and BAM formats for input and output.
- Configurable error rates and indel extension probabilities. Homopolymer length errors grow with run length.
- Optional chimeric reads.
- Fast: generates about 1 Gb of reads in less than 5 seconds on a 16 core CPU.

## Motivation

Oxford Nanopore data quality depends on many factors, such as the kit used, basecalling model version, and model precision level.
Basecalling models are updated often, so a simulator with fixed parameters goes out of date quickly.

ReadFaker takes read lengths and quality strings from your own data, so the simulated reads resembles your sequenced data, without need for tuning.
Errors are then added according to those quality scores.
The error model has defaults measured on recent data (see [Error Model Defaults](#error-model-defaults)), and all of them can be changed.

This is particularly useful for artificially contaminating real data for testing purposes (the reason I wrote this tool to begin with).

## Installation

Go to the [Releases](https://github.com/dialvarezs/readfaker/releases) and download the latest binary for Linux or Windows.
For other platforms, see [Building from Source](#building-from-source).

## Usage

```bash
readfaker -r <reference> -i <input> -o <output> -n <num_reads>
```

### Required Arguments

- `-r, --reference <FASTA>` - Reference sequences to sample reads from (optionally gzip-compressed)
- `-i, --input <FILE>` - Input file to extract quality and length models (FASTQ or BAM; alignments in the BAM are ignored)
- `-o, --output <FILE>` - Output file for simulated reads (FASTQ or unaligned BAM, detected by extension)

### Optional Arguments

- `-n, --num-reads <N>` - Number of reads to generate (default: 100000)
- `-s, --seed <N>` - Random seed for reproducibility
- `--compression-threads <N>` - Number of compression threads for output (default: 4, use 0 to detect them automatically)
- `--compression-level <N>` - BGZF compression level for output, 0-12 (default: 2)
- `--error-sub <RATE>` - Fraction of errors that are substitutions (default: 0.55)
- `--error-ins <RATE>` - Fraction of errors that are insertions (default: 0.18)
- `--error-del <RATE>` - Fraction of errors that are deletions (default: 0.27)
- `--error-ins-ext <RATE>` - Probability of extending an insertion by one more base (default: 0.27)
- `--error-del-ext <RATE>` - Probability of extending a deletion by one more base (default: 0.3)
- `--circular` - Treat all reference sequences as circular, so reads can span the end and start of a contig (e.g. bacterial chromosomes and plasmids)
- `--chimera-rate <RATE>` - Fraction of reads that are chimeras, joining two fragments from unrelated positions (default: 0, see [Chimeric Reads](#chimeric-reads))
- `--no-origin` - Don't write each read's origin (see [Read Origin](#read-origin))
- `-v, --verbose` - Enable verbose output

The three `--error-sub`, `--error-ins` and `--error-del` values set the mix of error types, not how many errors there are, and must add up to 1 or less.
How many errors a read gets depends on its quality scores.

### Examples

```bash
# Generate 10000 reads with verbose output
readfaker -r genome.fasta -i real_reads.fastq.gz -o simulated_reads.fastq.gz -n 10000 -v

# Generate reproducible reads with a fixed seed
readfaker -r genome.fasta -i real_reads.fastq -o simulated_reads.fastq -s 42

# Use BAM input and output with custom error rates
readfaker -r genome.fasta -i real_reads.bam -o simulated_reads.bam -n 50000 --error-sub 0.6 --error-ins 0.15 --error-del 0.25

# Adjust compression threads for better performance
readfaker -r genome.fasta -i real_reads.fastq.gz -o simulated_reads.fastq.gz -n 1000000 --compression-threads 8
```

## How It Works

1. Reads the input FASTQ or BAM file and stores its read lengths and quality strings.
2. Loads the reference sequences. Soft-masked (lowercase) bases are uppercased, and non-ACGT bases (N, IUPAC codes) are replaced by a random base.
3. For each read, samples a length, takes a fragment of that length from a random position and strand of the reference, and assigns it a real quality string of similar length.
4. Adds errors. Each base has an error with the probability given by its quality score, and the error is a substitution, insertion or deletion according to the configured mix. Homopolymer runs (3+ identical bases) also get length errors, mostly shortening, which are more likely in longer runs and in lower-quality reads.
5. Writes FASTQ or BAM records, with multi-threaded BGZF compression for `.gz`, `.bgz`, or `.bam` files.

## Error Model Defaults

Default error rates and the homopolymer model were measured on R10.4.1 reads (dorado v2.1 fast, hac and sup) from bacterial isolates aligned to their own assemblies, and HG002 aligned to its Q100 assembly.

## Read Origin

Reads are sampled from both strands with equal probability.
Each read records where it was sampled from as an `XO:Z:<contig>:<start>-<end>:<strand>` tag (1-based, inclusive reference coordinates; strand `+` or `-`).
With `--circular`, a read that spans the end of a contig has `end < start` (e.g. `plasmid1:4801-300:+`).
In FASTQ output it is the header comment; in BAM output it is an aux tag.
Use `--no-origin` to omit it, e.g. for blind tests where reads must not be told apart from real ones.

```
@5272f39c-f9a9-4f93-9cce-849905073826 XO:Z:chr1:3686-4184:-
```

## Chimeric Reads

With `--chimera-rate`, that fraction of reads joins two fragments taken from unrelated positions, contigs and strands.
The read length is split between both fragments, so read lengths still follow the input distribution.
The origin lists both fragments separated by `;`:

```
@f03d282b-a558-4bab-b35e-958fe339e183 XO:Z:chr1:2167214-2168091:+;chr2:1874138-1875122:-
```

The rate depends on the library prep, so chimeras are off by default.
In the data used for the [error model defaults](#error-model-defaults), about 3% of reads from a ligation library (HG002) and about 0.2% from a rapid barcoding library (bacterial isolates) were chimeras.

## Limitations

- Quality scores are taken at face value. Real basecallers are not perfectly calibrated, so simulated reads can have more or less substitutions than real ones.
- Junk reads, adapters and barcodes are not simulated.

## Building from Source

```bash
cargo build --release
```

The binary will be available at `target/release/readfaker`.

## Running Tests

```bash
cargo test
```