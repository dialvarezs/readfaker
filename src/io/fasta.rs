//! FASTA file reading using needletail.

use anyhow::{Context, Result, bail};
use needletail::{Sequence, parse_fastx_file};
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

        let mut reader = parse_fastx_file(path)
            .with_context(|| format!("Failed to open FASTA file: {}", path.display()))?;

        while let Some(record) = reader.next() {
            let record = record
                .with_context(|| format!("Failed to parse FASTA record in {}", path.display()))?;

            let id = String::from_utf8_lossy(record.id()).to_string();
            let sequence = record.normalize(false).to_vec();

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
