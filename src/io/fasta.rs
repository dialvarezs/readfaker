//! FASTA file reading.

use anyhow::{Context, Result, bail};
use noodles::fasta;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Represents a FASTA sequence with its ID and nucleotide sequence.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: Vec<u8>,
}

/// Reader for FASTA files.
pub struct FastaReader;

impl FastaReader {
    /// Reads all sequences from a FASTA file.
    ///
    /// # Arguments
    /// * `path` - Path to the FASTA file
    ///
    /// # Returns
    /// Vector of all FASTA records in the file
    pub fn read(path: &Path) -> Result<Vec<FastaRecord>> {
        let mut records = Vec::new();

        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTA file: {}", path.display()))?;
        let mut reader = fasta::io::Reader::new(BufReader::new(file));

        for result in reader.records() {
            let record = result
                .with_context(|| format!("Failed to parse FASTA record in {}", path.display()))?;

            let id = String::from_utf8_lossy(record.name()).to_string();
            let sequence = record.sequence().as_ref().to_vec();

            records.push(FastaRecord { id, sequence });
        }

        if records.is_empty() {
            bail!("No sequences found in FASTA file: {}", path.display());
        }

        Ok(records)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_record_creation() {
        let record = FastaRecord {
            id: "seq1".to_string(),
            sequence: b"ACGT".to_vec(),
        };
        assert_eq!(record.id, "seq1");
        assert_eq!(record.sequence, b"ACGT");
    }
}
