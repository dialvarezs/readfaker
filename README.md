# ReadFaker

A tool for simulating Oxford Nanopore sequencing reads with realistic quality profiles by extracting empirical models
from real FASTQ data.

## Features

- Creates empirical models for read length and quality scores (quality scores are grouped by length batches).
- Supports compressed input and output FASTQ files.
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

- Insertions and deletions are limited to one nucleotide length. Alteration type (substitution, insertion, deletion) ratios are fixed.
- Only generates modified sequences, not chimeras, junk reads and other types of artifacts.
- No BAM files support.

## Installation

Go to the [Releases](https://github.com/dialvarezs/readfaker/releases) and download the latest binary for your
platform.

## Usage

```bash
readfaker -r <reference.fasta> -i <input.fastq> -o <output.fastq> -n <num_reads>
```

### Required Arguments

- `-r, --reference <FASTA>` - Reference sequences to sample reads from
- `-i, --input <FASTQ>` - Input FASTQ file to extract quality and length models
- `-o, --output <FASTQ>` - Output FASTQ file for simulated reads

### Optional Arguments

- `-n, --num-reads <N>` - Number of reads to generate (default: 100000)
- `-s, --seed <N>` - Random seed for reproducibility
- `-v, --verbose` - Enable verbose output

### Example

```bash
# Generate 10000 reads with verbose output
readfaker -r genome.fasta -i real_reads.fastq.gz -o simulated_reads.fastq.gz -n 10000 -v

# Generate reproducible reads with a fixed seed
readfaker -r genome.fasta -i real_reads.fastq -o simulated_reads.fastq -s 42
```

## How It Works

1. **Model Extraction**: Reads an existing FASTQ file to build empirical models of read lengths and quality scores
2. **Reference Loading**: Parses reference genome sequences from FASTA format
3. **Read Generation**: Samples read lengths, selects random reference positions, applies quality profiles, and
   introduces errors based on quality scores
4. **Output**: Writes FASTQ records with automatic BGZF compression for `.gz`, `.bgz`, or `.bgzf` files

## Building from Source

```bash
cargo build --release
```

The binary will be available at `target/release/readfaker`.

## Running Tests

```bash
cargo test
```