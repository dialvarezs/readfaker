//! FASTQ file reading and writing.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use noodles::bgzf;
use noodles::fastq;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Reader for FASTQ files
pub struct FastqReader {
    reader: fastq::io::Reader<Box<dyn BufRead>>,
}

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
    ///     println!("Read: {}, length: {}", std::str::from_utf8(record.name()).unwrap(), record.sequence().len());
    /// }
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTQ file: {}", path.display()))?;

        // Check if file is gzip-compressed by reading magic bytes
        let mut buffered = BufReader::new(file);
        let is_compressed = is_gzip_compressed(&mut buffered)?;

        let reader: Box<dyn BufRead> = if is_compressed {
            // Use MultiGzDecoder which handles both regular gzip and BGZF
            Box::new(BufReader::new(MultiGzDecoder::new(buffered)))
        } else {
            Box::new(buffered)
        };

        Ok(Self {
            reader: fastq::io::Reader::new(reader),
        })
    }
}

impl Iterator for FastqReader {
    type Item = Result<fastq::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = fastq::Record::default();

        match self.reader.read_record(&mut record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(record)),
            Err(e) => Some(Err(
                anyhow::Error::new(e).context("Failed to parse FASTQ record")
            )),
        }
    }
}

/// Internal writer implementation supporting both uncompressed and BGZF-compressed output.
enum FastqWriterInner {
    Uncompressed(fastq::io::Writer<BufWriter<File>>),
    Compressed(fastq::io::Writer<bgzf::io::MultithreadedWriter<File>>),
}

impl FastqWriterInner {
    fn write_record(&mut self, record: &fastq::Record) -> std::io::Result<()> {
        match self {
            FastqWriterInner::Uncompressed(w) => w.write_record(record),
            FastqWriterInner::Compressed(w) => w.write_record(record),
        }
    }

    fn flush_writer(&mut self) -> std::io::Result<()> {
        match self {
            FastqWriterInner::Uncompressed(w) => w.get_mut().flush(),
            FastqWriterInner::Compressed(w) => w.get_mut().flush(),
        }
    }

    fn finish(self) -> Result<()> {
        match self {
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
/// use readfaker::io::fastq::FastqWriter;
/// use noodles::fastq;
/// use std::path::PathBuf;
///
/// // Uncompressed output
/// let mut writer = FastqWriter::new(&PathBuf::from("output.fastq"), 4)?;
///
/// // Compressed output (BGZF) with auto-detected threads
/// let mut writer_gz = FastqWriter::new(&PathBuf::from("output.fastq.gz"), 0)?;
///
/// // Compressed output with 4 threads
/// let mut writer_4t = FastqWriter::new(&PathBuf::from("output.fastq.gz"), 4)?;
///
/// let record = fastq::Record::new(
///     fastq::record::Definition::new("read1", ""),
///     b"ACGT",
///     b"IIII",
/// );
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
    /// * `compression_threads` - Number of compression threads (0 = auto-detect)
    pub fn new(path: &PathBuf, compression_threads: usize) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create FASTQ file: {}", path.display()))?;

        let writer = if should_compress(path) {
            // Use specified threads or auto-detect CPU cores
            let worker_count = if compression_threads == 0 {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(4) // Fallback to 4 threads
            } else {
                compression_threads
            };

            let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(
                std::num::NonZero::new(worker_count).unwrap(),
                file,
            );

            FastqWriterInner::Compressed(fastq::io::Writer::new(bgzf_writer))
        } else {
            FastqWriterInner::Uncompressed(fastq::io::Writer::new(BufWriter::new(file)))
        };

        Ok(Self { writer })
    }

    /// Writes a single FASTQ record.
    pub fn write_record(&mut self, record: &fastq::Record) -> Result<()> {
        self.writer
            .write_record(record)
            .context("Failed to write FASTQ record")
    }

    /// Writes multiple FASTQ records.
    ///
    /// # Arguments
    /// * `records` - Slice of FASTQ records to write
    pub fn write_records(&mut self, records: &[fastq::Record]) -> Result<()> {
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
        self.writer
            .flush_writer()
            .context("Failed to flush FASTQ writer")
    }

    /// Finishes the writer, properly shutting down compression threads if applicable.
    ///
    /// For compressed writers, this shuts down the thread pool and writes the final BGZF EOF block.
    /// This should be called explicitly before the writer is dropped to ensure proper finalization.
    pub fn finish(self) -> Result<()> {
        self.writer.finish()
    }
}

/// Helper function to check if a file is gzip-compressed
fn is_gzip_compressed<R: std::io::Read>(reader: &mut BufReader<R>) -> Result<bool> {
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
    fn test_fastq_writer() {
        let temp_dir = std::env::temp_dir();
        let temp_file = temp_dir.join("test.fastq");

        {
            let mut writer = FastqWriter::new(&temp_file, 4).unwrap();
            let record = fastq::Record::new(
                fastq::record::Definition::new("read1", ""),
                b"ACGT",
                b"IIII",
            );
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
            let mut writer = FastqWriter::new(&temp_file, 4).unwrap();
            let record1 = fastq::Record::new(
                fastq::record::Definition::new("read1", ""),
                b"ACGT",
                b"IIII",
            );
            let record2 = fastq::Record::new(
                fastq::record::Definition::new("read2", ""),
                b"TGCATGCA",
                b"IIIIIIII",
            );

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
        let records: Vec<fastq::Record> = FastqReader::from_path(&temp_file)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), "read1");
        assert_eq!(records[0].sequence(), b"ACGT");
        assert_eq!(records[1].name(), "read2");

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
