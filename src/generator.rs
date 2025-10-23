use crate::io::fasta::FastaRecord;
use crate::io::fastq::FastqReader;
use crate::io::FastqRecord;
use crate::models::{ErrorModel, LengthModel, QualityModel};
use crate::utils::{get_random_nucleotide, QUALITY_MAPPING};
use anyhow::{anyhow, Result};
use rand::prelude::IndexedRandom;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::path::Path;

const PHRED_OFFSET: u8 = 33;

/// Generator for synthetic sequencing reads with realistic error profiles.
///
/// Produces FASTQ records by sampling subsequences from reference genomes and applying
/// sequencing errors based on quality score models. Designed for simulating
/// contamination reads in genome assembly exercises.
///
/// # Example
/// ```no_run
/// use readfaker::generator::ReadGenerator;
/// use readfaker::models::{LengthModel, QualityModel};
/// use readfaker::io::fasta::FastaRecord;
///
/// let references = vec![FastaRecord {
///     id: "ecoli".to_string(),
///     sequence: b"ACGTACGT".to_vec(),
/// }];
/// let mut length_dist = LengthModel::new();
/// length_dist.add_value(100);
/// let mut quality_dist = QualityModel::new();
/// quality_dist.add_value(100, vec![b'I'; 100]);
///
/// let mut generator = ReadGenerator::new(references, length_dist, quality_dist, Some(42)).unwrap();
/// let read = generator.generate_read().unwrap();
/// ```
pub struct ReadGenerator {
    reference_sequences: Vec<FastaRecord>,
    length_model: LengthModel,
    quality_model: QualityModel,
    error_model: ErrorModel,
    rng: StdRng,
}

impl ReadGenerator {
    /// Creates a new read generator with specified models and random seed.
    ///
    /// # Arguments
    /// * `reference_sequences` - Reference genomes to sample subsequences from (must not be empty)
    /// * `length_model` - Empirical model of read lengths
    /// * `quality_model` - Empirical model of quality scores by read length
    /// * `seed` - Optional random seed for reproducibility (uses system entropy if None)
    ///
    /// # Returns
    /// A configured `ReadGenerator` ready to produce reads
    ///
    /// # Errors
    /// Returns an error if `reference_sequences` is empty
    pub fn new(
        reference_sequences: Vec<FastaRecord>,
        length_model: LengthModel,
        quality_model: QualityModel,
        error_model: ErrorModel,
        seed: Option<u64>,
    ) -> Result<Self> {
        if reference_sequences.is_empty() {
            anyhow::bail!("Reference sequences cannot be empty");
        }

        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_rng(&mut rand::rng()),
        };

        Ok(Self {
            reference_sequences,
            length_model,
            quality_model,
            error_model,
            rng,
        })
    }

    /// Generates a single synthetic read with realistic sequencing errors.
    ///
    /// Samples a read length from the model, chooses a random reference sequence,
    /// extracts a random subsequence, applies quality-based errors, and returns a FASTQ record.
    /// Automatically retries if the sampled length exceeds the reference sequence length.
    ///
    /// # Returns
    /// A `FastqRecord` with simulated sequencing errors based on quality scores
    ///
    /// # Errors
    /// Returns an error if the length or quality models are empty
    pub fn generate_read(&mut self) -> Result<FastqRecord> {
        loop {
            let length = self
                .length_model
                .sample(&mut self.rng)
                .ok_or_else(|| anyhow!("Length model is empty"))?;
            let reference_sequence = self.reference_sequences.choose(&mut self.rng).unwrap();

            // Skip if sampled length is longer than reference sequence
            if length > reference_sequence.sequence.len() {
                continue;
            }

            let max_start = reference_sequence.sequence.len() - length;
            let start_position = self.rng.random_range(0..=max_start);
            let mut sequence =
                reference_sequence.sequence[start_position..start_position + length].to_vec();

            let Some(qualities) = self.quality_model.sample(length, &mut self.rng) else {
                continue; // Skip if no quality string available
            };

            // Insert sequence nucleotide substitutions depending on error value
            qualities
                .iter()
                .enumerate()
                .for_each(|(i, &quality_ascii)| {
                    let phred = usize::from(quality_ascii.saturating_sub(PHRED_OFFSET).min(93));
                    let error_probability = QUALITY_MAPPING[phred];
                    if self.rng.random_range(0.0..1.0) <= error_probability {
                        sequence[i] = get_random_nucleotide(sequence[i], &mut self.rng);
                    }
                });

            return Ok(FastqRecord {
                id: format!("read_{}", self.rng.random::<u64>()),
                sequence,
                quality: qualities,
            });
        }
    }
}

/// Loads length and quality models from an existing FASTQ file.
///
/// Reads all records from the input file and builds empirical models
/// for read lengths and quality scores.
///
/// # Arguments
/// * `fastq_path` - Path to the FASTQ file to analyze
///
/// # Returns
/// Tuple of (LengthModel, QualityModel) built from the input file
pub fn load_models(fastq_path: &Path) -> Result<(LengthModel, QualityModel)> {
    let mut length_model = LengthModel::new();
    let mut quality_model = QualityModel::new();

    let reader = FastqReader::from_path(fastq_path)?;

    for record in reader {
        let record = record?;
        length_model.add_value(record.len());
        quality_model.add_value(record.len(), record.quality);
    }

    Ok((length_model, quality_model))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_read() {
        let sequences = vec![FastaRecord {
            id: "seq1".to_string(),
            sequence: b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        }];

        let mut length_model = LengthModel::new();
        let mut quality_model = QualityModel::new();
        let error_model = ErrorModel::new(None, None, None);
        length_model.add_value(10);
        quality_model.add_value(10, vec![b'?'; 10]); // Phred 30 as ASCII

        let mut generator = ReadGenerator::new(
            sequences,
            length_model,
            quality_model,
            error_model,
            Some(42),
        )
        .unwrap();

        // Generate multiple reads to verify the generator can be reused
        for _ in 0..5 {
            let read = generator.generate_read().unwrap();
            assert_eq!(read.len(), 10);
            assert!(read.quality.iter().all(|&q| q >= PHRED_OFFSET));
            assert!(read.id.starts_with("read_"));
        }
    }

    #[test]
    fn test_empty_reference_sequences() {
        let sequences = vec![]; // Empty!
        let length_model = LengthModel::new();
        let quality_model = QualityModel::new();
        let error_model = ErrorModel::new(None, None, None);

        let result = ReadGenerator::new(
            sequences,
            length_model,
            quality_model,
            error_model,
            Some(42),
        );
        assert!(result.is_err());
        let err = result.err().unwrap();
        assert_eq!(err.to_string(), "Reference sequences cannot be empty");
    }
}
