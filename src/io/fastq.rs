//! FASTQ file reading and writing.

use anyhow::{Context, Result};
use bgzip::Compression;
use bgzip::write::BGZFMultiThreadWriter;
use needletail::{Sequence, parse_fastx_file};
use std::fs::File;
use std::io::{self, BufWriter, Write};
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

    /// Checks whether the sequence is empty.
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

/// Internal writer implementation supporting both uncompressed and BGZF-compressed output.
enum FastqWriterInner {
    Uncompressed(BufWriter<File>),
    Compressed(BGZFMultiThreadWriter<File>),
}

impl Write for FastqWriterInner {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            FastqWriterInner::Uncompressed(w) => w.write(buf),
            FastqWriterInner::Compressed(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            FastqWriterInner::Uncompressed(w) => w.flush(),
            // Don't flush BGZFWriter - it handles finalization in its Drop implementation
            FastqWriterInner::Compressed(_) => Ok(()),
        }
    }
}

/// Writer for FASTQ files supporting both uncompressed and BGZF-compressed output.
///
/// Writes FASTQ records using buffered I/O. When using the `new()` constructor,
/// compression is automatically enabled based on the file extension (`.gz`, `.bgz`, or `.bgzf`).
///
/// The buffer is automatically flushed when the writer is dropped, but flush
/// errors are silently ignored. Call `flush()` explicitly if you need to
/// handle potential I/O errors during the final flush.
///
/// # Example
/// ```no_run
/// use readfaker::io::fastq::{FastqWriter, FastqRecord};
/// use std::path::PathBuf;
///
/// // Uncompressed output
/// let mut writer = FastqWriter::new(&PathBuf::from("output.fastq"))?;
///
/// // Compressed output (BGZF)
/// let mut writer_gz = FastqWriter::new(&PathBuf::from("output.fastq.gz"))?;
///
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ACGT".to_vec(),
///     b"IIII".to_vec()
/// )?;
/// writer.write_record(&record)?;
/// writer.flush()?;  // Explicitly flush to handle errors
/// # Ok::<(), anyhow::Error>(())
/// ```
pub struct FastqWriter {
    writer: FastqWriterInner,
}

impl FastqWriter {
    /// Creates a new FASTQ writer for the specified file path.
    ///
    /// The output format is automatically determined from the file extension:
    /// - Files ending in `.gz`, `.bgz`, or `.bgzf` will be BGZF-compressed
    /// - All other files will be uncompressed
    ///
    /// # Arguments
    /// * `path` - Path to the output FASTQ file
    pub fn new(path: &PathBuf) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create FASTQ file: {}", path.display()))?;

        let writer = if should_bgzf_compress(path) {
            FastqWriterInner::Compressed(BGZFMultiThreadWriter::new(file, Compression::default()))
        } else {
            FastqWriterInner::Uncompressed(BufWriter::new(file))
        };

        Ok(Self { writer })
    }

    /// Writes a single FASTQ record.
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        (|| -> io::Result<()> {
            writeln!(self.writer, "@{}", record.id)?;
            self.writer.write_all(&record.sequence)?;
            self.writer.write_all(b"\n+\n")?;
            self.writer.write_all(&record.quality)?;
            self.writer.write_all(b"\n")?;
            Ok(())
        })()
        .context("Failed to write FASTQ record")
    }

    /// Writes multiple FASTQ records.
    ///
    /// # Arguments
    /// * `records` - Slice of FASTQ records to write
    pub fn write_records(&mut self, records: &[FastqRecord]) -> Result<()> {
        for record in records {
            self.write_record(record)?;
        }
        Ok(())
    }

    /// Flushes the internal buffer to ensure all data is written.
    ///
    /// It's recommended to call this explicitly before the writer is dropped
    /// to ensure flush errors are properly handled.
    ///
    /// **Note:** For BGZF-compressed writers, this does not force immediate finalization.
    /// The compressed stream is properly closed when the writer is dropped.
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("Failed to flush FASTQ writer")
    }
}

fn should_bgzf_compress(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| {
            ["gz", "bgz", "bgzf"]
                .iter()
                .any(|s| ext.eq_ignore_ascii_case(s))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastq_record_new() {
        let record = FastqRecord::new("read1".to_string(), b"ACGT".to_vec(), b"IIII".to_vec());
        assert!(record.is_ok());
        assert_eq!(record.unwrap().len(), 4);
    }

    #[test]
    fn test_fastq_writer() {
        let temp_dir = std::env::temp_dir();
        let temp_file = temp_dir.join("test.fastq");

        {
            let mut writer = FastqWriter::new(&temp_file).unwrap();
            let record =
                FastqRecord::new("read1".to_string(), b"ACGT".to_vec(), b"IIII".to_vec()).unwrap();
            writer.write_record(&record).unwrap();
            writer.flush().unwrap();
        }

        let content = std::fs::read_to_string(&temp_file).unwrap();
        assert!(content.contains("@read1"));
        assert!(content.contains("ACGT"));

        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_fastq_writer_compressed() {
        let temp_file = std::env::temp_dir().join("readfaker_test_compressed.fastq.gz");
        std::fs::remove_file(&temp_file).ok();

        {
            let mut writer = FastqWriter::new(&temp_file).unwrap();
            let record1 =
                FastqRecord::new("read1".to_string(), b"ACGT".to_vec(), b"IIII".to_vec()).unwrap();
            let record2 = FastqRecord::new(
                "read2".to_string(),
                b"TGCATGCA".to_vec(),
                b"IIIIIIII".to_vec(),
            )
            .unwrap();

            writer.write_record(&record1).unwrap();
            writer.write_record(&record2).unwrap();
            writer.flush().unwrap();
        }

        // Verify the file is actually gzip-compressed by checking magic bytes
        let mut file = File::open(&temp_file).unwrap();
        let mut magic_bytes = [0u8; 2];
        io::Read::read_exact(&mut file, &mut magic_bytes).unwrap();
        assert_eq!(
            magic_bytes,
            [0x1f, 0x8b],
            "File should have gzip magic bytes"
        );

        // Read back and verify content
        let records: Vec<FastqRecord> = FastqReader::from_path(&temp_file)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "read1");
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[1].id, "read2");

        std::fs::remove_file(temp_file).ok();
    }
    #[test]
    fn test_should_bgzf_compress_suffixes() {
        assert!(should_bgzf_compress(Path::new("reads.fastq.gz")));
        assert!(should_bgzf_compress(Path::new("reads.fastq.bgz")));
        assert!(should_bgzf_compress(Path::new("reads.fastq.bgzf")));
        assert!(should_bgzf_compress(Path::new("reads.GZ")));
        assert!(!should_bgzf_compress(Path::new("reads.fastq")));
    }
}
