// I/O module for reading and writing sequence files

pub mod fasta;
pub mod fastq;

// Re-export main types
pub use fasta::FastaReader;
pub use fastq::{FastqRecord, FastqWriter};
