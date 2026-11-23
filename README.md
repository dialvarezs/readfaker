<p align="center"><img src="assets/logo.svg" alt="ReadFaker" width="360"></p>

# ReadFaker

A simple and fast Oxford Nanopore sequencing reads with realistic quality profiles by extracting empirical models from real sequencing data.

## Features

- Creates empirical models for read length and quality scores.
- Supports both FASTQ and BAM formats for input and output.
- Configurable error rates and indel extension probabilities. Homopolymer length errors grow with run length.
- Optional chimeric reads.
- Fast: generates about 1 Gb of reads in less than 5 seconds.

## Motivation

Oxford Nanopore data quality depends on many factors, such as the kit used, basecalling model version, and model precision level.
Basecalling models are updated often, so a simulator with fixed parameters goes out of date quickly.

ReadFaker takes read lengths and quality strings from your own data, so the simulated reads resembles your sequenced data, without need for tuning.
Errors are then added according to those quality scores.
The error model has defaults measured on recent data (see [Error Model Defaults](#error-model-defaults)), and all of them can be changed.

This is particularly useful for artificially contaminating real data for testing purposes (the reason I wrote this tool to begin with).

## Installation

ReadFaker is available on Bioconda:

```
conda install -c conda-forge -c bioconda readfaker
```

Pre-built binaries for Linux, Windows and macOS are on the [Releases](https://github.com/dialvarezs/readfaker/releases) page.
For other platforms, see [Building from Source](#building-from-source).

A container image is also available for use in pipelines, at [`ghcr.io/dialvarezs/readfaker`](https://github.com/dialvarezs/readfaker/pkgs/container/readfaker).

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
- `-a, --abundance <TSV>` - Relative abundance of each reference sequence (see [Abundance](#abundance))
- `-s, --seed <N>` - Random seed for reproducibility (read names included)
- `--compression-threads <N>` - Number of compression threads for output (default: 4, use 0 to detect them automatically)
- `--compression-level <N>` - BGZF compression level for output, 0-12 (default: 2)
- `--error-sub <RATE>` - Fraction of errors that are substitutions (default: 0.59)
- `--error-ins <RATE>` - Fraction of errors that are insertions (default: 0.18)
- `--error-del <RATE>` - Fraction of errors that are deletions (default: 0.23)
- `--error-ins-ext <RATE>` - Probability of extending an insertion by one more base (default: 0.26)
- `--error-del-ext <RATE>` - Probability of extending a deletion by one more base (default: 0.34)
- `--circular` - Treat all reference sequences as circular, so reads can span the end and start of a contig (e.g. bacterial chromosomes and plasmids)
- `--chimera-rate <RATE>` - Fraction of reads that are chimeras, joining two fragments from unrelated positions (default: 0, see [Chimeric Reads](#chimeric-reads))
- `--no-origin` - Don't write each read's origin (see [Read Origin](#read-origin))
- `-v, --verbose` - Enable verbose output

> [!NOTE]
> The three `--error-sub`, `--error-ins` and `--error-del` values set the mix of error types outside homopolymer runs, not how many errors there are, and must add up to 1.

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
4. Adds errors. The quality scores give the number of wrong bases expected in the read. Homopolymer runs (3+ identical bases) get length errors first, more likely in longer runs and in lower-quality reads. The remaining expected errors are spread over the other bases in proportion to their error probability, as substitutions, insertions or deletions according to the configured mix. Inside homopolymer runs only substitutions are added this way.
5. Writes FASTQ or BAM records, with multi-threaded BGZF compression for `.gz`, `.bgz`, or `.bam` files.

### Abundance

By default, reads are spread over the reference sequences in proportion to their length, so all of them get the same coverage.
With `--abundance`, each sequence gets a relative abundance (sequencing depth), and is picked in proportion to its length times its abundance.
This is useful to simulate a metagenome from several genomes.
The file has two tab-separated columns, the sequence id (as in the FASTA header, up to the first space) and its abundance.
Every reference sequence must be listed, and an abundance of 0 leaves a sequence out.
Lines starting with `#` are ignored.

```
# genome A at 10x the depth of genome B; A's plasmid follows its chromosome
genomeA_chr	10
genomeA_plasmid	10
genomeB_chr	1
```

### Read Origin

Reads are sampled from both strands with equal probability.
Each read records where it was sampled from as an `XO:Z:<contig>:<start>-<end>:<strand>` tag (1-based, inclusive reference coordinates; strand `+` or `-`).
With `--circular`, a read that spans the end of a contig has `end < start` (e.g. `plasmid1:4801-300:+`).
In FASTQ output it is the header comment; in BAM output it is an aux tag.
Use `--no-origin` to omit it, e.g. for blind tests where reads must not be told apart from real ones.

```
@5272f39c-f9a9-4f93-9cce-849905073826 XO:Z:chr1:3686-4184:-
```

### Chimeric Reads

With `--chimera-rate`, that fraction of reads joins two fragments taken from unrelated positions, contigs and strands.
The read length is split between both fragments, so read lengths still follow the input distribution.
The origin lists both fragments separated by `;`:

```
@f03d282b-a558-4bab-b35e-958fe339e183 XO:Z:chr1:2167214-2168091:+;chr2:1874138-1875122:-
```

The rate depends on the library prep, so chimeras are off by default.
In the data used for the [error model defaults](#error-model-defaults), about 3% of reads from a ligation library (HG002) and about 0.2% from a rapid barcoding library (bacterial isolates) were chimeras.

## Error Model Defaults

Default error rates and the homopolymer model were measured on R10.4.1 reads (dorado fast@v5.2.0, hac@v6.0.0 and sup@v5.2.0) from bacterial isolates aligned to their own assemblies, and HG002 aligned to its Q100 assembly.

## Limitations

- Quality scores are taken at face value, so simulated reads are only as accurate as the basecaller's quality scores. In some datasets, reads had between 11% fewer and 12% more errors than their quality scores predict. On older models this can be worse.
- Homopolymer errors depend on the base (G/C runs have about twice the errors of A/T runs) and vary between datasets at the same read quality, but the model only uses run length and read quality.
- Junk reads, adapters and barcodes are not simulated.

## Comparison with Other Simulators

ReadFaker needs no training step, reproduces homopolymer errors more closely than other simulators, and is much faster.

|                        | ReadFaker 0.4.0 | Badread 0.4.2 | NanoSim 3.2.3 | PBSIM3 3.0.5 (sample) |
|------------------------|-----------------|---------------|---------------|-----------------------|
| Training needed        | No              | Yes (33 min)  | Yes (5 min)   | No                    |
| Homopolymer error      | **0.018**       | 0.045         | 0.103         | 0.095                 |
| Read identity error    | **0.001**       | 0.005         | 0.003         | 0.007                 |
| Time per Gb (1 thread) | **3 s**         | 44 min        | 4 min         | 3 min                 |
| Peak memory            | 96 MB           | 417 MB        | 4.2 GB        | **33 MB**             |

- Homopolymer error: mean absolute difference from real reads in the fraction of homopolymer runs with a wrong length, over run lengths 3 to 9.
- Read identity error: absolute difference from real reads in mean read identity.

Accuracy was measured on R10.4.1 reads of the ZymoBIOMICS mock community (dorado fast, hac and sup), which were not used to set ReadFaker's defaults. Badread and NanoSim were trained on the same reads they were compared against. Time and memory are from a single run generating 1 Gb of human reads.

## Building from Source

```bash
cargo build --release
```

The binary will be available at `target/release/readfaker`.
