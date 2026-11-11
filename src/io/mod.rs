//! I/O module for reading and writing sequence files.
//!
//! Provides readers and writers for FASTA, FASTQ, and BAM file formats.

pub mod bam;
pub mod fasta;
pub mod fastq;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use noodles::bgzf;
pub use noodles::bgzf::io::writer::CompressionLevel;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::NonZero;
use std::path::Path;

// Re-export main types
pub use bam::{BamReader, BamWriter};
pub use fasta::FastaReader;
pub use fastq::FastqWriter;

/// Creates a multithreaded BGZF writer (`compression_threads` = 0 auto-detects CPU cores).
fn bgzf_writer(
    file: File,
    compression_threads: usize,
    compression_level: CompressionLevel,
) -> bgzf::io::MultithreadedWriter<File> {
    let worker_count = NonZero::new(compression_threads)
        .or_else(|| std::thread::available_parallelism().ok())
        .unwrap_or(NonZero::new(4).unwrap());
    bgzf::io::multithreaded_writer::Builder::default()
        .set_worker_count(worker_count)
        .set_compression_level(compression_level)
        .build_from_writer(file)
}

/// Opens a file for buffered reading, transparently decompressing gzip/BGZF (detected by magic bytes).
fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let mut buffered = BufReader::new(file);
    let header = buffered.fill_buf().context("Failed to read file header")?;

    // MultiGzDecoder handles both regular gzip and BGZF
    Ok(if header.starts_with(&[0x1f, 0x8b]) {
        Box::new(BufReader::new(MultiGzDecoder::new(buffered)))
    } else {
        Box::new(buffered)
    })
}
