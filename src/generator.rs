use crate::io::fasta::FastaRecord;
use crate::models::error::AlterationType;
use crate::models::{ErrorModel, LengthModel, QualityModel};
use crate::utils::{PHRED_OFFSET, QUALITY_MAPPING, make_rng};
use anyhow::{Result, anyhow, bail};
use noodles::fastq;
use rand::RngExt;
use rand::distr::Distribution;
use rand::distr::weighted::WeightedIndex;
use rand::rngs::StdRng;
use uuid::Uuid;

/// Maximum attempts to find a (length, reference, quality) combination before giving up.
const MAX_ATTEMPTS: usize = 1000;

/// Generator for synthetic sequencing reads with realistic error profiles.
///
/// Produces FASTQ records by sampling subsequences from reference genomes and applying
/// sequencing errors based on quality score models.
///
/// # Example
/// ```no_run
/// use readfaker::generator::ReadGenerator;
/// use readfaker::models::{ErrorModel, LengthModel, QualityModel};
/// use readfaker::io::fasta::FastaRecord;
/// use rand::SeedableRng;
/// use rand::rngs::StdRng;
///
/// let references = vec![FastaRecord {
///     id: "ecoli".to_string(),
///     sequence: b"ACGTACGT".to_vec(),
/// }];
/// let mut length_model = LengthModel::new();
/// length_model.add_value(100);
/// let mut quality_model = QualityModel::new(None, None, None);
/// let mut rng = StdRng::seed_from_u64(42);
/// quality_model.add_value(100, vec![b':'; 100], &mut rng);
/// let error_model = ErrorModel::new(None, None, None, None, None).unwrap();
///
/// let mut generator = ReadGenerator::new(
///     references,
///     length_model,
///     quality_model,
///     error_model,
///     Some(42),
/// ).unwrap();
/// let read = generator.generate_read().unwrap();
/// ```
pub struct ReadGenerator {
    reference_sequences: Vec<FastaRecord>,
    /// Picks references proportionally to their length, so coverage is uniform across contigs.
    reference_weights: WeightedIndex<usize>,
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
            bail!("Reference sequences cannot be empty");
        }

        let reference_weights =
            WeightedIndex::new(reference_sequences.iter().map(|r| r.sequence.len()))
                .map_err(|e| anyhow!("Invalid reference sequences: {}", e))?;
        let rng = make_rng(seed);

        Ok(Self {
            reference_sequences,
            reference_weights,
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
    /// A `fastq::Record` with simulated sequencing errors based on quality scores
    ///
    /// # Errors
    /// Returns an error if the length model is empty, or if no valid read can be built
    /// after `MAX_ATTEMPTS` tries (e.g. references shorter than all sampled lengths)
    pub fn generate_read(&mut self) -> Result<fastq::Record> {
        for _ in 0..MAX_ATTEMPTS {
            let length = self
                .length_model
                .sample(&mut self.rng)
                .ok_or_else(|| anyhow!("Length model is empty"))?;
            let reference_sequence =
                &self.reference_sequences[self.reference_weights.sample(&mut self.rng)];

            // Skip if sampled length is longer than reference sequence
            if length > reference_sequence.sequence.len() {
                continue;
            }

            let max_start = reference_sequence.sequence.len() - length;
            let start_position = self.rng.random_range(0..=max_start);
            let sequence =
                reference_sequence.sequence[start_position..start_position + length].to_vec();

            let Some(qualities) = self.quality_model.sample(length, &mut self.rng) else {
                continue; // Skip if no quality string available
            };

            let (final_sequence, final_qualities) = self.apply_errors(sequence, qualities);

            let id = Uuid::new_v4().to_string();
            return Ok(fastq::Record::new(
                fastq::record::Definition::new(id, ""),
                final_sequence,
                final_qualities,
            ));
        }

        bail!(
            "Failed to generate a read after {} attempts: sampled read lengths may exceed all reference sequences",
            MAX_ATTEMPTS
        )
    }

    /// Applies sequencing errors to a sequence based on quality scores and error model.
    ///
    /// For each position, uses the quality score to determine if an error occurs,
    /// then uses the error model to determine the type of error (substitution, insertion, or deletion).
    ///
    /// # Arguments
    /// * `sequence` - Original nucleotide sequence
    /// * `qualities` - Quality scores for each position
    ///
    /// # Returns
    /// Tuple of (modified sequence, modified quality scores)
    fn apply_errors(&mut self, sequence: Vec<u8>, qualities: Vec<u8>) -> (Vec<u8>, Vec<u8>) {
        let mut new_sequence = Vec::with_capacity(sequence.len());
        let mut new_quality = Vec::with_capacity(qualities.len());

        let mut i = 0;
        while i < sequence.len() {
            let quality_ascii = qualities[i];
            let phred = usize::from(quality_ascii.saturating_sub(PHRED_OFFSET).min(93));
            let error_probability = QUALITY_MAPPING[phred];
            let alteration = if self.rng.random_range(0.0..1.0) <= error_probability {
                self.error_model.get_alteration_type(&mut self.rng)
            } else {
                None
            };

            match alteration {
                Some(AlterationType::Substitution) => {
                    new_sequence.push(self.get_random_nucleotide(Some(sequence[i])));
                    new_quality.push(quality_ascii);
                }
                Some(AlterationType::Insertion(count)) => {
                    new_sequence.push(sequence[i]);
                    new_quality.push(quality_ascii);

                    for _ in 0..count {
                        new_sequence.push(self.get_random_nucleotide(None));
                        new_quality.push(quality_ascii); // reuse the same quality for inserted bases
                    }
                }
                Some(AlterationType::Deletion(count)) => {
                    let skip = count.max(1).min(sequence.len() - i);
                    i += skip;
                    continue;
                }
                None => {
                    new_sequence.push(sequence[i]);
                    new_quality.push(quality_ascii);
                }
            }

            i += 1;
        }

        (new_sequence, new_quality)
    }

    /// Returns a random nucleotide, optionally excluding a specific one.
    ///
    /// # Arguments
    /// * `exclude` - Optional nucleotide to exclude (as ASCII byte: b'A', b'C', b'G', or b'T')
    ///   - `Some(nucleotide)` returns a different nucleotide
    ///   - `None` returns any random nucleotide
    ///
    /// # Returns
    /// A random nucleotide byte from {A, C, G, T}
    fn get_random_nucleotide(&mut self, exclude: Option<u8>) -> u8 {
        const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];

        match exclude {
            Some(nucleotide) => {
                let idx = NUCLEOTIDES
                    .iter()
                    .position(|&n| n == nucleotide)
                    .unwrap_or(0);
                let offset = self.rng.random_range(1..=3);
                NUCLEOTIDES[(idx + offset) % 4]
            }
            None => {
                let idx = self.rng.random_range(0..4);
                NUCLEOTIDES[idx]
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn create_test_generator(sequences: Option<Vec<FastaRecord>>) -> Result<ReadGenerator> {
        let sequences = sequences.unwrap_or_else(|| {
            vec![FastaRecord {
                id: "seq1".to_string(),
                sequence: b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            }]
        });

        let mut length_model = LengthModel::new();
        let mut quality_model = QualityModel::new(None, None, None);
        let error_model = ErrorModel::new(None, None, None, None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        length_model.add_value(10);
        quality_model.add_value(10, vec![b'?'; 10], &mut rng); // Phred 30 as ASCII

        ReadGenerator::new(
            sequences,
            length_model,
            quality_model,
            error_model,
            Some(42),
        )
    }

    #[test]
    fn test_empty_reference_sequences() {
        let result = create_test_generator(Some(vec![]));
        assert!(result.is_err());
        let err = result.err().unwrap();
        assert_eq!(err.to_string(), "Reference sequences cannot be empty");
    }

    #[test]
    fn test_generate_read() {
        let mut generator = create_test_generator(None).unwrap();

        // Generate multiple reads to verify the generator can be reused
        for _ in 0..5 {
            let read = generator.generate_read().unwrap();
            assert_eq!(read.sequence().len(), 10);
            assert!(read.quality_scores().iter().all(|&q| q >= PHRED_OFFSET));
            // Verify the ID is a valid UUID
            let name_str = std::str::from_utf8(read.name()).unwrap();
            assert!(
                Uuid::parse_str(name_str).is_ok(),
                "Expected valid UUID, got: {}",
                name_str
            );
        }
    }

    #[test]
    fn test_get_random_nucleotide() {
        let mut generator = create_test_generator(None).unwrap();

        // Test with exclusion (substitution)
        let result = generator.get_random_nucleotide(Some(b'A'));
        assert_ne!(result, b'A');
        assert!(result == b'C' || result == b'G' || result == b'T');

        // Test without exclusion (insertion)
        let result = generator.get_random_nucleotide(None);
        assert!(result == b'A' || result == b'C' || result == b'G' || result == b'T');
    }

    #[test]
    fn test_generate_read_fails_when_reference_too_short() {
        let mut generator = create_test_generator(Some(vec![FastaRecord {
            id: "short".to_string(),
            sequence: b"ACGT".to_vec(),
        }]))
        .unwrap();
        assert!(generator.generate_read().is_err());
    }

    #[test]
    fn test_reference_choice_weighted_by_length() {
        let long = FastaRecord {
            id: "long".to_string(),
            sequence: vec![b'A'; 990],
        };
        let short = FastaRecord {
            id: "short".to_string(),
            sequence: vec![b'C'; 10],
        };
        let mut generator = create_test_generator(Some(vec![long, short])).unwrap();

        let from_short = (0..1000)
            .filter(|_| generator.generate_read().unwrap().sequence()[0] == b'C')
            .count();
        // Expected ~1% (10/1000 bp); uniform choice would give ~50%
        assert!(
            from_short < 50,
            "got {} reads from short contig",
            from_short
        );
    }
}
