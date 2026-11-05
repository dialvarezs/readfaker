//! I/O module for reading and writing sequence files.
//!
//! Provides readers and writers for FASTA, FASTQ, and BAM file formats.

pub mod bam;
pub mod fasta;
pub mod fastq;

use noodles::bgzf;
use std::fs::File;
use std::num::NonZero;

// Re-export main types
pub use bam::{BamReader, BamWriter};
pub use fasta::FastaReader;
pub use fastq::FastqWriter;

/// Creates a multithreaded BGZF writer (`compression_threads` = 0 auto-detects CPU cores).
fn bgzf_writer(file: File, compression_threads: usize) -> bgzf::io::MultithreadedWriter<File> {
    let worker_count = NonZero::new(compression_threads)
        .or_else(|| std::thread::available_parallelism().ok())
        .unwrap_or(NonZero::new(4).unwrap());
    bgzf::io::MultithreadedWriter::with_worker_count(worker_count, file)
}
