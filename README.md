# ReadFaker

A tool for simulating Oxford Nanopore sequencing reads with realistic quality profiles by extracting empirical models
from real FASTQ data.

## Features

- Creates empirical models for read length and quality scores (quality scores are grouped by length batches).
- Supports both FASTQ and BAM formats for input and output.
- Automatic compression detection for input files (gzip/BGZF).
- Multi-threaded BGZF compression for output files.
- Configurable error rates and indel extension probabilities.
- Fast: can generate a million reads in under a minute.

## Motivation

Oxford Nanopore data quality depends on many factors, such as the kit used, basecalling model version, and model
precision level.
Basecalling models keep improving quite often, making it challenging to simulate realistic data with fixed parameters.

This tool takes a different approach: instead of hardcoded models, it extracts length and quality profiles directly from
your real data, ensuring the simulated reads match the characteristics of actual sequencing runs.

This is particularly useful for artificially contaminating real data for testing purposes (the reason I wrote this tool
to begin with).

## Current Limitations / Planned Improvements

- Only generates modified sequences, not chimeras, junk reads and other types of artifacts.

## Installation

Go to the [Releases](https://github.com/dialvarezs/readfaker/releases) and download the latest binary for your
platform.

## Usage

```bash
readfaker -r <reference> -i <input> -o <output> -n <num_reads>
```

### Required Arguments

- `-r, --reference <FASTA>` - Reference sequences to sample reads from
- `-i, --input <FILE>` - Input file to extract quality and length models (FASTQ or BAM)
- `-o, --output <FILE>` - Output file for simulated reads (FASTQ or BAM, detected by extension)

### Optional Arguments

- `-n, --num-reads <N>` - Number of reads to generate (default: 100000)
- `-s, --seed <N>` - Random seed for reproducibility
- `--compression-threads <N>` - Number of compression threads for output (default: 4)
- `--error-sub <RATE>` - Error substitution rate (default: 0.7)
- `--error-ins <RATE>` - Error insertion rate (default: 0.1)
- `--error-del <RATE>` - Error deletion rate (default: 0.2)
- `--error-ins-ext <RATE>` - Insertion extension probability using geometric distribution (default: 0.4)
- `--error-del-ext <RATE>` - Deletion extension probability using geometric distribution (default: 0.4)
- `-v, --verbose` - Enable verbose output

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

1. **Model Extraction**: Reads an existing FASTQ or BAM file to build empirical models of read lengths and quality scores
2. **Reference Loading**: Parses reference genome sequences from FASTA format
3. **Read Generation**: Samples read lengths, selects random reference positions, applies quality profiles, and
   introduces errors based on quality scores with configurable error rates and indel extension probabilities
4. **Output**: Writes FASTQ or BAM records with automatic multi-threaded BGZF compression for `.gz`, `.bgz`, `.bgzf`, or `.bam` files

## Building from Source

```bash
cargo build --release
```

The binary will be available at `target/release/readfaker`.

## Running Tests

```bash
cargo test
```