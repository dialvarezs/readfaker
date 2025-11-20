//! FASTQ file reading and writing.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use noodles::bgzf;
use noodles::fastq;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
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
    /// Automatically detects and handles both compressed (gzip/bgzf) and uncompressed FASTQ files.
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
    pub fn from_path(path: &Path) -> Result<FastqReaderIterator> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTQ file: {}", path.display()))?;

        // Check if file is gzip-compressed by reading magic bytes
        let mut buffered = BufReader::new(file);
        let is_compressed = is_gzip_compressed(&mut buffered)?;

        let reader = if is_compressed {
            // Use MultiGzDecoder which handles both regular gzip and BGZF
            FastqReaderInner::Compressed(fastq::io::Reader::new(BufReader::new(
                MultiGzDecoder::new(buffered),
            )))
        } else {
            FastqReaderInner::Uncompressed(fastq::io::Reader::new(buffered))
        };

        Ok(FastqReaderIterator { reader })
    }
}

/// Internal reader implementation supporting both compressed and uncompressed files
enum FastqReaderInner<R: std::io::Read> {
    Uncompressed(fastq::io::Reader<BufReader<R>>),
    Compressed(fastq::io::Reader<BufReader<MultiGzDecoder<BufReader<R>>>>),
}

/// Iterator over FASTQ records
pub struct FastqReaderIterator {
    reader: FastqReaderInner<File>,
}

impl Iterator for FastqReaderIterator {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = fastq::Record::default();

        let read_result = match &mut self.reader {
            FastqReaderInner::Uncompressed(reader) => reader.read_record(&mut record),
            FastqReaderInner::Compressed(reader) => reader.read_record(&mut record),
        };

        match read_result {
            Ok(0) => None,
            Ok(_) => {
                let id = record.name().to_string();
                let sequence = record.sequence().to_vec();
                let quality = record.quality_scores().to_vec();

                Some(Ok(FastqRecord {
                    id,
                    sequence,
                    quality,
                }))
            }
            Err(e) => Some(Err(
                anyhow::Error::new(e).context("Failed to parse FASTQ record")
            )),
        }
    }
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
    /// * `compression_threads` - Optional number of compression threads (None = auto-detect)
    pub fn new(path: &PathBuf, compression_threads: Option<usize>) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create FASTQ file: {}", path.display()))?;

        let writer = if should_compress(path) {
            // Use specified threads or auto-detect CPU cores
            let worker_count = compression_threads.unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(4) // Fallback to 4 threads
            });

            let bgzf_writer = if let Some(count) = std::num::NonZero::new(worker_count) {
                bgzf::io::MultithreadedWriter::with_worker_count(count, file)
            } else {
                bgzf::io::MultithreadedWriter::new(file)
            };

            FastqWriterInner::Compressed(fastq::io::Writer::new(bgzf_writer))
        } else {
            FastqWriterInner::Uncompressed(fastq::io::Writer::new(BufWriter::new(file)))
        };

        Ok(Self { writer })
    }

    /// Writes a single FASTQ record.
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        let noodles_record = fastq::Record::new(
            fastq::record::Definition::new(record.id.clone(), ""),
            record.sequence.clone(),
            record.quality.clone(),
        );

        match &mut self.writer {
            FastqWriterInner::Uncompressed(w) => w.write_record(&noodles_record),
            FastqWriterInner::Compressed(w) => w.write_record(&noodles_record),
        }
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
    pub fn flush(&mut self) -> Result<()> {
        match &mut self.writer {
            FastqWriterInner::Uncompressed(w) => w.get_mut().flush(),
            FastqWriterInner::Compressed(w) => w.get_mut().flush(),
        }
        .context("Failed to flush FASTQ writer")
    }

    /// Finishes the writer, properly shutting down compression threads if applicable.
    ///
    /// For compressed writers, this shuts down the thread pool and writes the final BGZF EOF block.
    /// This should be called explicitly before the writer is dropped to ensure proper finalization.
    pub fn finish(self) -> Result<()> {
        match self.writer {
            FastqWriterInner::Uncompressed(mut w) => w
                .get_mut()
                .flush()
                .context("Failed to flush uncompressed writer"),
            FastqWriterInner::Compressed(w) => {
                // Get the underlying BGZF writer and finish it to shutdown threads and write EOF
                w.into_inner()
                    .finish()
                    .map(|_| ()) // Discard the returned File handle
                    .map_err(|e| anyhow::anyhow!("Failed to finish BGZF writer: {}", e))
            }
        }
    }
}

/// Internal writer implementation supporting both uncompressed and BGZF-compressed output.
enum FastqWriterInner {
    Uncompressed(fastq::io::Writer<BufWriter<File>>),
    Compressed(fastq::io::Writer<bgzf::io::MultithreadedWriter<File>>),
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
/// let mut writer = FastqWriter::new(&PathBuf::from("output.fastq"), None)?;
///
/// // Compressed output (BGZF) with auto-detected threads
/// let mut writer_gz = FastqWriter::new(&PathBuf::from("output.fastq.gz"), None)?;
///
/// // Compressed output with 4 threads
/// let mut writer_4t = FastqWriter::new(&PathBuf::from("output.fastq.gz"), Some(4))?;
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

/// Helper function to check if a file is gzip-compressed
fn is_gzip_compressed<R: std::io::Read>(reader: &mut BufReader<R>) -> Result<bool> {
    use std::io::BufRead;

    let buffer = reader.fill_buf().context("Failed to read file header")?;

    // Check for gzip magic bytes (0x1f 0x8b)
    Ok(buffer.len() >= 2 && buffer[0] == 0x1f && buffer[1] == 0x8b)
}

/// Helper function to check if a file should be compressed based on its extension
fn should_compress(path: &Path) -> bool {
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
            let mut writer = FastqWriter::new(&temp_file, None).unwrap();
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
            let mut writer = FastqWriter::new(&temp_file, None).unwrap();
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
        use std::io::Read;
        let mut file = File::open(&temp_file).unwrap();
        let mut magic_bytes = [0u8; 2];
        file.read_exact(&mut magic_bytes).unwrap();
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
        assert!(should_compress(Path::new("reads.fastq.gz")));
        assert!(should_compress(Path::new("reads.fastq.bgz")));
        assert!(should_compress(Path::new("reads.fastq.bgzf")));
        assert!(should_compress(Path::new("reads.GZ")));
        assert!(!should_compress(Path::new("reads.fastq")));
    }
}
