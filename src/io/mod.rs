//! I/O module for reading and writing sequence files.
//!
//! Provides readers and writers for FASTA, FASTQ, and BAM file formats.

pub mod bam;
pub mod fasta;
pub mod fastq;

// Re-export main types
pub use bam::{BamReader, BamWriter};
pub use fasta::FastaReader;
pub use fastq::FastqWriter;
