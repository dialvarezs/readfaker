//! FASTA file reading using needletail.

use anyhow::{Context, Result};
use needletail::{parse_fastx_file, Sequence};
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
            anyhow::bail!("No sequences found in FASTA file: {}", path.display());
        }

        Ok(records)
    }

    /// Calculates the total length of all sequences.
    ///
    /// # Arguments
    /// * `records` - Slice of FASTA records
    ///
    /// # Returns
    /// Sum of all sequence lengths
    pub fn total_length(records: &[FastaRecord]) -> usize {
        records.iter().map(|r| r.sequence.len()).sum()
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

    #[test]
    fn test_total_length() {
        let records = vec![
            FastaRecord {
                id: "seq1".to_string(),
                sequence: b"ACGT".to_vec(),
            },
            FastaRecord {
                id: "seq2".to_string(),
                sequence: b"TGCA".to_vec(),
            },
        ];
        assert_eq!(FastaReader::total_length(&records), 8);
    }
}
