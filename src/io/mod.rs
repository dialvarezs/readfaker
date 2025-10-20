//! I/O module for reading and writing sequence files.
//!
//! Provides readers and writers for FASTA and FASTQ file formats.

pub mod fasta;
pub mod fastq;

// Re-export main types
pub use fasta::FastaReader;
pub use fastq::{FastqRecord, FastqWriter};
