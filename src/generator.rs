use crate::io::fasta::FastaRecord;
use crate::models::error::AlterationType;
use crate::models::{ErrorModel, LengthModel, QualityModel};
use crate::utils::{ERROR_THRESHOLDS, make_rng};
use anyhow::{Result, anyhow, bail};
use noodles::fastq;
use rand::RngExt;
use rand::distr::Distribution;
use rand::distr::weighted::WeightedIndex;
use rand::rngs::SmallRng;
use std::ops::Range;
use uuid::Uuid;

/// Maximum attempts to find a (length, reference, quality) combination before giving up.
const MAX_ATTEMPTS: usize = 1000;

/// SAM tag holding the read's origin (`contig:start-end:strand`, 1-based inclusive).
/// Written as the FASTQ description (`XO:Z:...`) and as a BAM aux tag.
pub const ORIGIN_TAG: &str = "XO";

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
    rng: SmallRng,
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
        mut reference_sequences: Vec<FastaRecord>,
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
        let mut rng = make_rng(seed);

        // Nanopore reads are uppercase ACGT only: uppercase soft-masked bases and
        // replace N/IUPAC codes with a random base.
        for record in &mut reference_sequences {
            for base in &mut record.sequence {
                *base = match base.to_ascii_uppercase() {
                    b @ (b'A' | b'C' | b'G' | b'T') => b,
                    _ => b"ACGT"[rng.random_range(0..4)],
                };
            }
        }

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
            let mut sequence =
                reference_sequence.sequence[start_position..start_position + length].to_vec();
            let reverse = self.rng.random_bool(0.5);
            if reverse {
                reverse_complement(&mut sequence);
            }
            let origin = format!(
                "{ORIGIN_TAG}:Z:{}:{}-{}:{}",
                reference_sequence.id,
                start_position + 1,
                start_position + length,
                if reverse { '-' } else { '+' }
            );

            let Some(qualities) = self.quality_model.sample(length, &mut self.rng) else {
                continue; // Skip if no quality string available
            };

            let (sequence, qualities) = self.apply_homopolymer_errors(sequence, qualities);
            let (final_sequence, final_qualities) = self.apply_errors(sequence, qualities);

            let id = Uuid::new_v4().to_string();
            return Ok(fastq::Record::new(
                fastq::record::Definition::new(id, origin),
                final_sequence,
                final_qualities,
            ));
        }

        bail!(
            "Failed to generate a read after {} attempts: sampled read lengths may exceed all reference sequences",
            MAX_ATTEMPTS
        )
    }

    /// Applies homopolymer length errors on top of the per-base errors.
    ///
    /// Each run of 3+ identical bases may be shortened or lengthened. The probability grows with
    /// the run length and with the read's mean error probability.
    ///
    /// # Arguments
    /// * `sequence` - Original nucleotide sequence
    /// * `qualities` - Quality scores for each position
    ///
    /// # Returns
    /// Tuple of (modified sequence, modified quality scores); the inputs are returned unchanged
    /// when no run is altered
    fn apply_homopolymer_errors(
        &mut self,
        sequence: Vec<u8>,
        qualities: Vec<u8>,
    ) -> (Vec<u8>, Vec<u8>) {
        // Use the read's mean: qualities come from a different read, so the values at the run's
        // positions say nothing about the run.
        let probabilities = self
            .error_model
            .homopolymer_error_probabilities(mean_error_probability(&qualities));

        // Bases before `copied` are already in the new vectors; nothing is copied until a run changes.
        let mut new_sequence = Vec::new();
        let mut new_quality = Vec::new();
        let mut copied = 0;

        for run in homopolymer_runs(&sequence) {
            let run_length = run.len();
            let probability = probabilities[run_length.min(probabilities.len() - 1)];
            if self.rng.random::<f64>() >= probability {
                continue;
            }

            new_sequence.extend_from_slice(&sequence[copied..run.end]);
            new_quality.extend_from_slice(&qualities[copied..run.end]);
            copied = run.end;

            let change = self
                .error_model
                .get_homopolymer_alteration(run_length, &mut self.rng);
            let count = change.unsigned_abs();
            if change < 0 {
                // keep at least one base of the run
                let keep = new_sequence.len() - count.min(run_length - 1);
                new_sequence.truncate(keep);
                new_quality.truncate(keep);
            } else {
                new_sequence.extend(std::iter::repeat_n(sequence[run.start], count));
                new_quality.extend(std::iter::repeat_n(qualities[run.end - 1], count));
            }
        }

        if copied == 0 {
            return (sequence, qualities);
        }
        new_sequence.extend_from_slice(&sequence[copied..]);
        new_quality.extend_from_slice(&qualities[copied..]);
        (new_sequence, new_quality)
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

        let thresholds = &*ERROR_THRESHOLDS;
        // Error-free stretches are copied in bulk; only error positions are handled per base
        let mut run_start = 0;
        let mut i = 0;
        while i < sequence.len() {
            let quality_ascii = qualities[i];
            let alteration = if self.rng.random::<u32>() <= thresholds[usize::from(quality_ascii)] {
                self.error_model.get_alteration_type(&mut self.rng)
            } else {
                None
            };
            let Some(alteration) = alteration else {
                i += 1;
                continue;
            };

            new_sequence.extend_from_slice(&sequence[run_start..i]);
            new_quality.extend_from_slice(&qualities[run_start..i]);

            match alteration {
                AlterationType::Substitution => {
                    new_sequence.push(self.get_random_nucleotide(Some(sequence[i])));
                    new_quality.push(quality_ascii);
                    i += 1;
                }
                AlterationType::Insertion(count) => {
                    new_sequence.push(sequence[i]);
                    new_quality.push(quality_ascii);

                    for _ in 0..count {
                        new_sequence.push(self.get_random_nucleotide(None));
                        new_quality.push(quality_ascii); // reuse the same quality for inserted bases
                    }
                    i += 1;
                }
                AlterationType::Deletion(count) => {
                    i += count.max(1).min(sequence.len() - i);
                }
            }
            run_start = i;
        }
        new_sequence.extend_from_slice(&sequence[run_start..]);
        new_quality.extend_from_slice(&qualities[run_start..]);

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

/// Reverse-complements an uppercase ACGT sequence in place.
///
/// # Arguments
/// * `sequence` - Nucleotide sequence to reverse-complement
fn reverse_complement(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = match *base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        };
    }
}

/// Finds all homopolymer runs of 3 or more identical bases.
///
/// Scans 64-base blocks with a branch-free bitmask (bit j set when bases j..j+3 are equal), since
/// a base-by-base run scan mispredicts at nearly every run boundary.
///
/// # Arguments
/// * `sequence` - Nucleotide sequence to scan
///
/// # Returns
/// Position ranges of the runs, in order
fn homopolymer_runs(sequence: &[u8]) -> Vec<Range<usize>> {
    let n = sequence.len();
    let mut runs: Vec<Range<usize>> = Vec::new();
    let mut block = 0;
    while block + 2 < n {
        let width = (n - 2 - block).min(64);
        let mut mask = 0u64;
        for j in 0..width {
            let k = block + j;
            let equal = (sequence[k] == sequence[k + 1]) & (sequence[k + 1] == sequence[k + 2]);
            mask |= u64::from(equal) << j;
        }

        while mask != 0 {
            let start = block + mask.trailing_zeros() as usize;
            mask &= mask - 1;
            if runs.last().is_some_and(|run| start < run.end) {
                continue; // inside the previous run
            }
            let mut end = start + 3;
            while end < n && sequence[end] == sequence[start] {
                end += 1;
            }
            runs.push(start..end);
        }
        block += 64;
    }
    runs
}

/// Computes the mean per-base error probability of a read from its quality scores.
///
/// # Arguments
/// * `qualities` - Quality scores (Phred+33 ASCII)
///
/// # Returns
/// Mean error probability, between 0 and 1
fn mean_error_probability(qualities: &[u8]) -> f64 {
    let thresholds = &*ERROR_THRESHOLDS;
    let sum: u64 = qualities
        .iter()
        .map(|&q| u64::from(thresholds[usize::from(q)]))
        .sum();
    sum as f64 / (qualities.len() as f64 * f64::from(u32::MAX))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

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
            assert!(
                read.quality_scores()
                    .iter()
                    .all(|&q| q >= crate::utils::PHRED_OFFSET)
            );
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
            .filter(|_| {
                generator
                    .generate_read()
                    .unwrap()
                    .description()
                    .starts_with(b"XO:Z:short:")
            })
            .count();
        // Expected ~1% (10/1000 bp); uniform choice would give ~50%
        assert!(
            from_short < 50,
            "got {} reads from short contig",
            from_short
        );
    }

    /// Generator whose quality strings are all Q93, so reads are effectively error-free.
    fn error_free_generator(sequence: &[u8], length: usize) -> ReadGenerator {
        let mut length_model = LengthModel::new();
        let mut quality_model = QualityModel::new(None, None, None);
        length_model.add_value(length);
        quality_model.add_value(length, vec![b'~'; length], &mut StdRng::seed_from_u64(1));
        ReadGenerator::new(
            vec![FastaRecord {
                id: "chr1".to_string(),
                sequence: sequence.to_vec(),
            }],
            length_model,
            quality_model,
            ErrorModel::new(None, None, None, None, None).unwrap(),
            Some(7),
        )
        .unwrap()
    }

    #[test]
    fn test_origin_matches_sequence_on_both_strands() {
        let reference = b"AACCGTTAGCATGCATTTGACCAGTAGGCATCAGGTACCGATCGATGCA";
        let mut generator = error_free_generator(reference, 12);
        let (mut forward, mut reverse) = (0, 0);

        for _ in 0..200 {
            let read = generator.generate_read().unwrap();
            let origin = std::str::from_utf8(read.description()).unwrap();
            let (range, strand) = origin
                .strip_prefix("XO:Z:chr1:")
                .unwrap()
                .split_once(':')
                .unwrap();
            let (start, end) = range.split_once('-').unwrap();
            let (start, end): (usize, usize) = (start.parse().unwrap(), end.parse().unwrap());

            let mut expected = reference[start - 1..end].to_vec();
            match strand {
                "+" => forward += 1,
                "-" => {
                    reverse += 1;
                    reverse_complement(&mut expected);
                }
                _ => panic!("bad strand in {origin}"),
            }
            assert_eq!(read.sequence(), expected.as_slice(), "origin {origin}");
        }
        assert!(
            forward > 50 && reverse > 50,
            "{forward} fwd / {reverse} rev"
        );
    }

    #[test]
    fn test_homopolymer_errors_mostly_shorten_long_runs() {
        let mut generator = create_test_generator(None).unwrap();
        // Q0 => error probability 1, so every run with a non-zero weight is altered
        let (mut shorter, mut longer) = (0, 0);
        for _ in 0..1000 {
            let (seq, qual) =
                generator.apply_homopolymer_errors(b"CAAAAAAG".to_vec(), vec![b'!'; 8]);
            assert_eq!(seq.len(), qual.len());
            assert!(seq.starts_with(b"CA") && seq.ends_with(b"AG"));
            assert!(seq[1..seq.len() - 1].iter().all(|&b| b == b'A'));
            match seq.len().cmp(&8) {
                std::cmp::Ordering::Less => shorter += 1,
                std::cmp::Ordering::Greater => longer += 1,
                std::cmp::Ordering::Equal => panic!("run of 6 left unchanged at Q0"),
            }
        }
        assert!(
            (700..800).contains(&shorter),
            "{shorter} shorter / {longer} longer"
        );

        // Runs of 1-2 bases are left to the per-base error model
        let (seq, _) = generator.apply_homopolymer_errors(b"ACCGTTA".to_vec(), vec![b'!'; 7]);
        assert_eq!(seq, b"ACCGTTA");
    }

    #[test]
    fn test_homopolymer_scan_finds_runs_at_block_edges() {
        fn rle(seq: &[u8]) -> Vec<(u8, usize)> {
            seq.chunk_by(|a, b| a == b)
                .map(|r| (r[0], r.len()))
                .collect()
        }
        // Runs at the start, straddling the 64-base block boundary, and at the very end
        let mut seq = b"GGG".to_vec();
        seq.extend(b"ACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTAC"); // to 62
        seq.extend(b"TTTT"); // 62..66
        seq.extend(b"ACGACGACGACGACGACG");
        seq.extend(b"CCC");
        let mut generator = create_test_generator(None).unwrap();

        let (out, _) = generator.apply_homopolymer_errors(seq.clone(), vec![b'!'; seq.len()]);
        let (before, after) = (rle(&seq), rle(&out));
        assert_eq!(before.len(), after.len());
        let changed: Vec<_> = before
            .iter()
            .zip(&after)
            .filter(|(b, a)| b != a)
            .map(|(b, _)| *b)
            .collect();
        assert_eq!(changed, vec![(b'G', 3), (b'T', 4), (b'C', 3)]);
    }

    #[test]
    fn test_reference_normalized_to_uppercase_acgt() {
        let mut generator = error_free_generator(b"acgtNNRYacgt", 12);
        let read = generator.generate_read().unwrap();
        assert!(read.sequence().iter().all(|b| b"ACGT".contains(b)));
    }
}
