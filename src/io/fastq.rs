// FASTQ file reading and writing

use anyhow::{Context, Result};
use needletail::{parse_fastx_file, Sequence};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Represents a FASTQ record
#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub id: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

impl FastqRecord {
    /// Create a new FASTQ record
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

    /// Get the length of the read
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if the record is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
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
    pub fn from_path(
        path: &Path,
    ) -> Result<impl Iterator<Item = Result<FastqRecord>> + '_> {
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
/// use std::path::Path;
///
/// let mut writer = FastqWriter::new(Path::new("output.fastq"))?;
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
    /// Create a new FASTQ writer for the specified file path
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create FASTQ file: {}", path.display()))?;

        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    /// Write a single FASTQ record to the file
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        writeln!(self.writer, "@{}", record.id)
            .context("Failed to write FASTQ header")?;

        self.writer.write_all(&record.sequence)
            .context("Failed to write FASTQ sequence")?;
        writeln!(self.writer).context("Failed to write newline after sequence")?;

        writeln!(self.writer, "+")
            .context("Failed to write FASTQ separator")?;

        self.writer.write_all(&record.quality)
            .context("Failed to write FASTQ quality scores")?;
        writeln!(self.writer).context("Failed to write newline after quality")?;

        Ok(())
    }

    /// Flush the internal buffer to ensure all data is written to disk
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()
            .context("Failed to flush FASTQ writer")
    }
}

impl Drop for FastqWriter {
    fn drop(&mut self) {
        // Best effort flush on drop - ignore errors since we can't handle them in Drop
        let _ = self.writer.flush();
    }
}
