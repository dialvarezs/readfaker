//! FASTQ file reading and writing.

use anyhow::{Context, Result};
use needletail::{parse_fastx_file, Sequence};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Represents a FASTQ record with sequence and quality scores.
#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub id: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

impl FastqRecord {
    /// Creates a new FASTQ record with validation.
    ///
    /// # Arguments
    /// * `id` - Read identifier
    /// * `sequence` - Nucleotide sequence
    /// * `quality` - Quality scores (must match sequence length)
    ///
    /// # Returns
    /// Validated FASTQ record or error if lengths don't match
    pub fn new(id: String, sequence: Vec<u8>, quality: Vec<u8>) -> Result<Self> {
        if sequence.len() != quality.len() {
            anyhow::bail!(
                "Sequence length ({}) does not match quality length ({})",
                sequence.len(),
                quality.len()
            );
        }
        Ok(Self {
            id,
            sequence,
            quality,
        })
    }

    /// Returns the length of the sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}

/// Reader for FASTQ files
pub struct FastqReader;

impl FastqReader {
    /// Open a FASTQ file and return a streaming iterator over records
    ///
    /// # Example
    /// ```no_run
    /// use readfaker::io::fastq::FastqReader;
    /// use std::path::Path;
    ///
    /// let reader = FastqReader::from_path(Path::new("input.fastq"))?;
    /// for record in reader {
    ///     let record = record?;
    ///     println!("Read: {}, length: {}", record.id, record.len());
    /// }
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn from_path(path: &Path) -> Result<impl Iterator<Item = Result<FastqRecord>> + '_> {
        let mut reader = parse_fastx_file(path)
            .with_context(|| format!("Failed to open FASTQ file: {}", path.display()))?;

        Ok(std::iter::from_fn(move || {
            reader.next().map(|result| {
                result
                    .context("Failed to parse FASTQ record")
                    .and_then(|record| {
                        let id = String::from_utf8_lossy(record.id()).to_string();
                        let sequence = record.normalize(false).to_vec();
                        let quality = record
                            .qual()
                            .context("Missing quality scores in FASTQ record")?
                            .to_vec();

                        Ok(FastqRecord {
                            id,
                            sequence,
                            quality,
                        })
                    })
            })
        }))
    }
}

/// Writer for FASTQ files
///
/// Writes FASTQ records to a file using buffered I/O.
/// The buffer is automatically flushed when the writer is dropped.
///
/// # Example
/// ```no_run
/// use readfaker::io::fastq::{FastqWriter, FastqRecord};
/// use std::path::PathBuf;
///
/// let mut writer = FastqWriter::new(&PathBuf::from("output.fastq"))?;
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ACGT".to_vec(),
///     b"IIII".to_vec()
/// )?;
/// writer.write_record(&record)?;
/// writer.flush()?;
/// # Ok::<(), anyhow::Error>(())
/// ```
pub struct FastqWriter {
    writer: BufWriter<File>,
}

impl FastqWriter {
    /// Creates a new FASTQ writer for the specified file path.
    ///
    /// # Arguments
    /// * `path` - Path to the output FASTQ file
    pub fn new(path: &PathBuf) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create FASTQ file: {}", path.display()))?;

        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    /// Write a single FASTQ record to the file
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        (|| -> std::io::Result<()> {
            writeln!(self.writer, "@{}", record.id)?;
            self.writer.write_all(&record.sequence)?;
            self.writer.write_all(b"\n+\n")?;
            self.writer.write_all(&record.quality)?;
            self.writer.write_all(b"\n")?;
            Ok(())
        })()
        .context("Failed to write FASTQ record")
    }

    /// Writes multiple FASTQ records to the file.
    ///
    /// # Arguments
    /// * `records` - Slice of FASTQ records to write
    pub fn write_records(&mut self, records: &[FastqRecord]) -> Result<()> {
        for record in records {
            self.write_record(record)?;
        }
        Ok(())
    }

    /// Flush the internal buffer to ensure all data is written to disk
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("Failed to flush FASTQ writer")
    }
}

impl Drop for FastqWriter {
    fn drop(&mut self) {
        let _ = self.writer.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastq_record_new() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ACGT".to_vec(),
            b"IIII".to_vec(),
        );
        assert!(record.is_ok());
        assert_eq!(record.unwrap().len(), 4);
    }

    #[test]
    fn test_fastq_writer() {
        let temp_dir = std::env::temp_dir();
        let temp_file = temp_dir.join("test.fastq");

        {
            let mut writer = FastqWriter::new(&temp_file).unwrap();
            let record = FastqRecord::new(
                "read1".to_string(),
                b"ACGT".to_vec(),
                b"IIII".to_vec(),
            ).unwrap();
            writer.write_record(&record).unwrap();
            writer.flush().unwrap();
        }

        let content = std::fs::read_to_string(&temp_file).unwrap();
        assert!(content.contains("@read1"));
        assert!(content.contains("ACGT"));

        std::fs::remove_file(temp_file).ok();
    }
}
