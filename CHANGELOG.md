# Changelog

## 0.3.0 - 2026-11-10

### Improvements

- Sample reads from both strands.
- Record the origin of each read (contig, coordinates and strand) as an `XO` tag in FASTQ and BAM output, disablable with `--no-origin`.
- Add homopolymer length errors, more likely in longer runs and in lower-quality reads.
- Add `--circular` to treat reference sequences as circular, so reads can span the end and start of a contig.
- Add `--chimera-rate` to generate chimeric reads that join two fragments from unrelated positions.
- Add `--compression-level` and use libdeflate for BGZF output.
- Support gzip-compressed FASTA references.
- Faster read generation (per-base error loop and read length sampling).

### Fixes

- Choose reference sequences in proportion to their length, so coverage is uniform across contigs.
- Stop with an error when no read can be generated, e.g. when all sampled lengths exceed the reference sequences.
- Cap the length of insertions and deletions.

### Changes

- New default error rates, based on R10.4.1 reads basecalled with dorado v2.1.2 fast, hac and sup.
- Reference sequences are uppercased, and non-ACGT bases are replaced by a random base.
- Output differs from 0.2.0 for the same seed.

## 0.2.0 - 2025-12-09

### Improvements

- Add support for BAM input/output.
- Implement multi-nucleotide insertions and deletions using a geometric distribution.
- Expose sequencing error rates as configurable CLI arguments.
- Introduce a CLI option to set the number of compression threads.

### Changes

- Migrated from needletail/bgzip to the noodles crate for unified file handling.

## 0.1.0 - 2025-10-26

Initial release.
