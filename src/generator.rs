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
use std::collections::HashMap;
use std::ops::Range;

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
///     false,
///     0.0,
///     Some(42),
/// ).unwrap();
/// let read = generator.generate_read().unwrap();
/// ```
pub struct ReadGenerator {
    reference_sequences: Vec<FastaRecord>,
    /// Picks references proportionally to their length (times their abundance, if given),
    /// so coverage is uniform across contigs with the same abundance.
    reference_weights: WeightedIndex<f64>,
    length_model: LengthModel,
    quality_model: QualityModel,
    error_model: ErrorModel,
    /// Reads may span the end and start of each reference sequence.
    circular: bool,
    /// Fraction of reads built from two unrelated fragments.
    chimera_rate: f64,
    rng: SmallRng,
}

impl ReadGenerator {
    /// Creates a new read generator with specified models and random seed.
    ///
    /// # Arguments
    /// * `reference_sequences` - Reference genomes to sample subsequences from (must not be empty)
    /// * `length_model` - Empirical model of read lengths
    /// * `quality_model` - Empirical model of quality scores by read length
    /// * `circular` - Treat all reference sequences as circular
    /// * `chimera_rate` - Fraction of reads that join two fragments from unrelated positions
    /// * `seed` - Optional random seed for reproducibility (uses system entropy if None)
    ///
    /// # Returns
    /// A configured `ReadGenerator` ready to produce reads
    ///
    /// # Errors
    /// Returns an error if `reference_sequences` is empty or `chimera_rate` is not between 0 and 1
    pub fn new(
        mut reference_sequences: Vec<FastaRecord>,
        length_model: LengthModel,
        quality_model: QualityModel,
        error_model: ErrorModel,
        circular: bool,
        chimera_rate: f64,
        seed: Option<u64>,
    ) -> Result<Self> {
        if !(0.0..=1.0).contains(&chimera_rate) {
            bail!("Chimera rate must be between 0 and 1, got {chimera_rate}");
        }
        if reference_sequences.is_empty() {
            bail!("Reference sequences cannot be empty");
        }

        let reference_weights =
            WeightedIndex::new(reference_sequences.iter().map(|r| r.sequence.len() as f64))
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
            circular,
            chimera_rate,
            rng,
        })
    }

    /// Sets the relative abundance (sequencing depth) of each reference sequence.
    ///
    /// Sequences are then picked in proportion to their length times their abundance.
    ///
    /// # Errors
    /// Returns an error if a reference sequence has no abundance, an abundance names an unknown
    /// sequence, or all weights are zero
    pub fn with_abundances(mut self, abundances: &HashMap<String, f64>) -> Result<Self> {
        let weights = self
            .reference_sequences
            .iter()
            .map(|r| {
                abundances
                    .get(&r.id)
                    .map(|a| r.sequence.len() as f64 * a)
                    .ok_or_else(|| anyhow!("No abundance given for reference sequence '{}'", r.id))
            })
            .collect::<Result<Vec<_>>>()?;
        if let Some(unknown) = abundances
            .keys()
            .find(|id| !self.reference_sequences.iter().any(|r| &r.id == *id))
        {
            bail!("Abundance given for unknown reference sequence '{unknown}'");
        }
        self.reference_weights =
            WeightedIndex::new(weights).map_err(|e| anyhow!("Invalid abundances: {}", e))?;
        Ok(self)
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
            // A chimeric read splits its length between two fragments, so read lengths still
            // follow the length model.
            let chimeric =
                self.chimera_rate > 0.0 && length > 1 && self.rng.random_bool(self.chimera_rate);
            let first_length = if chimeric {
                self.rng.random_range(1..length)
            } else {
                length
            };

            // Skip if a fragment is longer than its reference sequence
            let Some((mut sequence, mut origin)) = self.sample_fragment(first_length) else {
                continue;
            };
            if chimeric {
                let Some((fragment, fragment_origin)) = self.sample_fragment(length - first_length)
                else {
                    continue;
                };
                sequence.extend(fragment);
                origin.push(';');
                origin.push_str(&fragment_origin);
            }
            let origin = format!("{ORIGIN_TAG}:Z:{origin}");

            let Some(qualities) = self.quality_model.sample(length, &mut self.rng) else {
                continue; // Skip if no quality string available
            };

            let errors = expected_errors(&qualities);
            let (sequence, qualities, runs, homopolymer_errors) =
                self.apply_homopolymer_errors(sequence, qualities, errors);
            let (final_sequence, final_qualities) =
                self.apply_errors(sequence, qualities, &runs, errors, homopolymer_errors);

            // Drawn from the seeded RNG so read names are reproducible too
            let id = uuid::Builder::from_random_bytes(self.rng.random())
                .into_uuid()
                .to_string();
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

    /// Extracts a fragment from a random position and strand of a random reference sequence.
    ///
    /// # Arguments
    /// * `length` - Length of the fragment
    ///
    /// # Returns
    /// * `Some((sequence, origin))` - Fragment and its origin as `<contig>:<start>-<end>:<strand>`
    /// * `None` - If the chosen reference sequence is shorter than `length`
    fn sample_fragment(&mut self, length: usize) -> Option<(Vec<u8>, String)> {
        let reference_sequence =
            &self.reference_sequences[self.reference_weights.sample(&mut self.rng)];
        let contig = &reference_sequence.sequence;
        if length > contig.len() {
            return None;
        }

        let max_start = if self.circular {
            contig.len() - 1
        } else {
            contig.len() - length
        };
        let start_position = self.rng.random_range(0..=max_start);
        let end_position = start_position + length;
        let mut sequence = if end_position <= contig.len() {
            contig[start_position..end_position].to_vec()
        } else {
            [
                &contig[start_position..],
                &contig[..end_position - contig.len()],
            ]
            .concat()
        };
        let reverse = self.rng.random_bool(0.5);
        if reverse {
            reverse_complement(&mut sequence);
        }
        let origin = format!(
            "{}:{}-{}:{}",
            reference_sequence.id,
            start_position + 1,
            (end_position - 1) % contig.len() + 1,
            if reverse { '-' } else { '+' }
        );
        Some((sequence, origin))
    }

    /// Applies homopolymer length errors, before the per-base errors.
    ///
    /// Each run of 3+ identical bases may be shortened or lengthened. The probability grows with
    /// the run length and with the read's mean error probability.
    ///
    /// # Arguments
    /// * `sequence` - Original nucleotide sequence
    /// * `qualities` - Quality scores for each position
    /// * `errors` - Expected number of wrong bases in the read, from `expected_errors`
    ///
    /// # Returns
    /// Tuple of (modified sequence, modified quality scores, positions of the runs in the modified
    /// sequence, expected number of bases changed by run length errors); the inputs are returned
    /// unchanged when no run is altered
    fn apply_homopolymer_errors(
        &mut self,
        sequence: Vec<u8>,
        qualities: Vec<u8>,
        errors: f64,
    ) -> (Vec<u8>, Vec<u8>, Vec<Range<usize>>, f64) {
        // Use the read's mean: qualities come from a different read, so the values at the run's
        // positions say nothing about the run.
        let probabilities = self
            .error_model
            .homopolymer_error_probabilities(errors / qualities.len().max(1) as f64);

        // Bases before `copied` are already in the new vectors; nothing is copied until a run changes.
        let mut new_sequence = Vec::new();
        let mut new_quality = Vec::new();
        let mut copied = 0;
        let runs = homopolymer_runs(&sequence);
        let mut new_runs = Vec::with_capacity(runs.len());
        let mut expected_errors = 0.0;
        // Bases added (positive) or removed (negative) so far, to place runs in the new sequence
        let mut offset = 0isize;

        for run in runs {
            let run_length = run.len();
            let probability = probabilities[run_length.min(probabilities.len() - 1)];
            expected_errors +=
                probability * self.error_model.mean_homopolymer_error_bases(run_length);
            let start = run.start.saturating_add_signed(offset);
            if self.rng.random::<f64>() >= probability {
                new_runs.push(start..start + run_length);
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
                let removed = count.min(run_length - 1);
                let keep = new_sequence.len() - removed;
                new_sequence.truncate(keep);
                new_quality.truncate(keep);
                new_runs.push(start..start + run_length - removed);
                offset -= removed as isize;
            } else {
                new_sequence.extend(std::iter::repeat_n(sequence[run.start], count));
                new_quality.extend(std::iter::repeat_n(qualities[run.end - 1], count));
                new_runs.push(start..start + run_length + count);
                offset += count as isize;
            }
        }

        if copied == 0 {
            return (sequence, qualities, new_runs, expected_errors);
        }
        new_sequence.extend_from_slice(&sequence[copied..]);
        new_quality.extend_from_slice(&qualities[copied..]);
        (new_sequence, new_quality, new_runs, expected_errors)
    }

    /// Applies sequencing errors to a sequence based on quality scores and error model.
    ///
    /// The quality scores give the read's expected number of wrong bases. Homopolymer length
    /// errors are part of that number, so per-base error probabilities are scaled to fill the rest,
    /// counting each indel by its length. For each position, the scaled quality score determines if an error occurs, and the
    /// error model its type (substitution, insertion, or deletion). Inside homopolymer runs only
    /// substitutions apply, since their length errors were already drawn.
    ///
    /// # Arguments
    /// * `sequence` - Nucleotide sequence
    /// * `qualities` - Quality scores for each position
    /// * `runs` - Positions of homopolymer runs in `sequence`, in order
    /// * `errors` - Expected number of wrong bases in the read, from `expected_errors`
    /// * `homopolymer_errors` - Expected number of bases changed by homopolymer length errors
    ///
    /// # Returns
    /// Tuple of (modified sequence, modified quality scores)
    fn apply_errors(
        &mut self,
        sequence: Vec<u8>,
        qualities: Vec<u8>,
        runs: &[Range<usize>],
        errors: f64,
        homopolymer_errors: f64,
    ) -> (Vec<u8>, Vec<u8>) {
        let mut new_sequence = Vec::with_capacity(sequence.len());
        let mut new_quality = Vec::with_capacity(qualities.len());

        let thresholds = self.scaled_thresholds(&qualities, runs, errors, homopolymer_errors);
        // Error-free stretches are copied in bulk; only error positions are handled per base
        let mut run_start = 0;
        let mut i = 0;
        while i < sequence.len() {
            let quality_ascii = qualities[i];
            if self.rng.random::<u32>() > thresholds[usize::from(quality_ascii)] {
                i += 1;
                continue;
            }
            let alteration = self.error_model.get_alteration_type(&mut self.rng);
            let run = runs.partition_point(|run| run.end <= i);
            let in_run = runs.get(run).is_some_and(|run| run.start <= i);
            if in_run && !matches!(alteration, AlterationType::Substitution) {
                i += 1;
                continue;
            }

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

    /// Per-base error thresholds for one read, scaled so that per-base and homopolymer errors
    /// together change as many bases as its quality scores predict.
    ///
    /// # Arguments
    /// * `qualities` - Quality scores of the read
    /// * `runs` - Positions of homopolymer runs, where only substitutions apply
    /// * `errors` - Expected number of wrong bases in the read, from `expected_errors`
    /// * `homopolymer_errors` - Expected number of bases changed by homopolymer length errors
    ///
    /// # Returns
    /// Thresholds indexed by quality ASCII value, as in `ERROR_THRESHOLDS`
    fn scaled_thresholds(
        &self,
        qualities: &[u8],
        runs: &[Range<usize>],
        total: f64,
        homopolymer_errors: f64,
    ) -> [u32; 256] {
        let in_runs: f64 = runs
            .iter()
            .map(|run| expected_errors(&qualities[run.clone()]))
            .sum();

        let model = &self.error_model;
        let target = total - homopolymer_errors;
        let unscaled =
            model.mean_error_bases() * (total - in_runs) + model.substitution_rate * in_runs;
        let scale = if unscaled > 0.0 {
            (target / unscaled).max(0.0)
        } else {
            0.0
        };
        ERROR_THRESHOLDS.map(|t| (f64::from(t) * scale).min(f64::from(u32::MAX)) as u32)
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
/// Scans 64-base blocks with a branch-free bitmask (bit j set when bases j..j+3 are equal).
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

/// Computes the expected number of wrong bases from quality scores (sum of error probabilities).
///
/// # Arguments
/// * `qualities` - Quality scores (Phred+33 ASCII)
///
/// # Returns
/// Expected number of wrong bases
fn expected_errors(qualities: &[u8]) -> f64 {
    let thresholds = &*ERROR_THRESHOLDS;
    let sum: u64 = qualities
        .iter()
        .map(|&q| u64::from(thresholds[usize::from(q)]))
        .sum();
    sum as f64 / f64::from(u32::MAX)
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
            false,
            0.0,
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
                uuid::Uuid::parse_str(name_str).is_ok(),
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

    #[test]
    fn test_reference_choice_weighted_by_abundance() {
        let a = FastaRecord {
            id: "a".to_string(),
            sequence: vec![b'A'; 500],
        };
        let b = FastaRecord {
            id: "b".to_string(),
            sequence: vec![b'C'; 500],
        };
        let abundances = HashMap::from([("a".to_string(), 9.0), ("b".to_string(), 1.0)]);
        let mut generator = create_test_generator(Some(vec![a, b]))
            .unwrap()
            .with_abundances(&abundances)
            .unwrap();

        let from_b = (0..1000)
            .filter(|_| {
                generator
                    .generate_read()
                    .unwrap()
                    .description()
                    .starts_with(b"XO:Z:b:")
            })
            .count();
        // Expected ~10% (equal lengths, 1/(9+1) abundance)
        assert!((50..150).contains(&from_b), "got {} reads from b", from_b);
    }

    #[test]
    fn test_abundances_must_match_references() {
        let missing = HashMap::from([("other".to_string(), 1.0)]);
        assert!(
            create_test_generator(None)
                .unwrap()
                .with_abundances(&missing)
                .is_err()
        );
        let extra = HashMap::from([("seq1".to_string(), 1.0), ("other".to_string(), 1.0)]);
        assert!(
            create_test_generator(None)
                .unwrap()
                .with_abundances(&extra)
                .is_err()
        );
        let zero = HashMap::from([("seq1".to_string(), 0.0)]);
        assert!(
            create_test_generator(None)
                .unwrap()
                .with_abundances(&zero)
                .is_err()
        );
    }

    #[test]
    fn test_same_seed_gives_same_reads() {
        let mut a = create_test_generator(None).unwrap();
        let mut b = create_test_generator(None).unwrap();
        for _ in 0..10 {
            assert_eq!(a.generate_read().unwrap(), b.generate_read().unwrap());
        }
    }

    /// Generator whose quality strings are all Q93, so reads are effectively error-free.
    fn error_free_generator(
        sequence: &[u8],
        length: usize,
        circular: bool,
        chimera_rate: f64,
    ) -> ReadGenerator {
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
            circular,
            chimera_rate,
            Some(7),
        )
        .unwrap()
    }

    #[test]
    fn test_origin_matches_sequence_on_both_strands() {
        let reference = b"AACCGTTAGCATGCATTTGACCAGTAGGCATCAGGTACCGATCGATGCA";
        let mut generator = error_free_generator(reference, 12, false, 0.0);
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
    fn test_circular_reads_span_contig_end() {
        let reference = b"AACCGTTAGCATGCATTTGA";
        let mut generator = error_free_generator(reference, 12, true, 0.0);
        let mut wrapped = 0;

        for _ in 0..300 {
            let read = generator.generate_read().unwrap();
            let origin = std::str::from_utf8(read.description()).unwrap();
            let (range, strand) = origin
                .strip_prefix("XO:Z:chr1:")
                .unwrap()
                .split_once(':')
                .unwrap();
            let (start, end) = range.split_once('-').unwrap();
            let (start, end): (usize, usize) = (start.parse().unwrap(), end.parse().unwrap());

            let mut expected = if start <= end {
                reference[start - 1..end].to_vec()
            } else {
                wrapped += 1;
                [&reference[start - 1..], &reference[..end]].concat()
            };
            if strand == "-" {
                reverse_complement(&mut expected);
            }
            assert_eq!(read.sequence(), expected.as_slice(), "origin {origin}");
        }
        // 11 of 20 start positions wrap around
        assert!(wrapped > 100, "{wrapped} reads wrapped");
    }

    #[test]
    fn test_chimeric_reads_join_two_fragments() {
        let reference = b"AACCGTTAGCATGCATTTGACCAGTAGGCATCAGGTACCGATCGATGCA";
        let mut generator = error_free_generator(reference, 12, false, 0.5);
        let mut chimeras = 0;

        for _ in 0..200 {
            let read = generator.generate_read().unwrap();
            let origin = std::str::from_utf8(read.description()).unwrap();
            let segments: Vec<&str> = origin.strip_prefix("XO:Z:").unwrap().split(';').collect();
            assert!(segments.len() <= 2, "origin {origin}");
            chimeras += usize::from(segments.len() == 2);

            let mut expected = Vec::new();
            for segment in segments {
                let (range, strand) = segment
                    .strip_prefix("chr1:")
                    .unwrap()
                    .split_once(':')
                    .unwrap();
                let (start, end) = range.split_once('-').unwrap();
                let (start, end): (usize, usize) = (start.parse().unwrap(), end.parse().unwrap());
                let mut fragment = reference[start - 1..end].to_vec();
                if strand == "-" {
                    reverse_complement(&mut fragment);
                }
                expected.extend(fragment);
            }
            assert_eq!(read.sequence().len(), 12, "origin {origin}");
            assert_eq!(read.sequence(), expected.as_slice(), "origin {origin}");
        }
        assert!((70..130).contains(&chimeras), "{chimeras} chimeras");
    }

    #[test]
    fn test_invalid_chimera_rate() {
        let result = ReadGenerator::new(
            vec![FastaRecord {
                id: "chr1".to_string(),
                sequence: b"ACGT".to_vec(),
            }],
            LengthModel::new(),
            QualityModel::new(None, None, None),
            ErrorModel::new(None, None, None, None, None).unwrap(),
            false,
            1.5,
            None,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_homopolymer_errors_mostly_shorten_long_runs() {
        let mut generator = create_test_generator(None).unwrap();
        // Q0 => error probability 1, so most runs are altered
        let (mut shorter, mut longer) = (0, 0);
        for _ in 0..2000 {
            let (seq, qual, runs, _) =
                generator.apply_homopolymer_errors(b"CAAAAAAG".to_vec(), vec![b'!'; 8], 8.0);
            assert_eq!(seq.len(), qual.len());
            assert!(seq.starts_with(b"CA") && seq.ends_with(b"AG"));
            assert!(seq[1..seq.len() - 1].iter().all(|&b| b == b'A'));
            assert_eq!(runs, vec![1..seq.len() - 1]);
            match seq.len().cmp(&8) {
                std::cmp::Ordering::Less => shorter += 1,
                std::cmp::Ordering::Greater => longer += 1,
                std::cmp::Ordering::Equal => {}
            }
        }
        let share = f64::from(shorter) / f64::from(shorter + longer);
        assert!(
            (0.67..0.77).contains(&share),
            "{shorter} shorter / {longer} longer"
        );

        // Runs of 1-2 bases are left to the per-base error model
        let (seq, _, runs, expected) =
            generator.apply_homopolymer_errors(b"ACCGTTA".to_vec(), vec![b'!'; 7], 7.0);
        assert_eq!(seq, b"ACCGTTA");
        assert!(runs.is_empty());
        assert_eq!(expected, 0.0);
    }

    #[test]
    fn test_homopolymer_errors_come_out_of_the_quality_budget() {
        let generator = create_test_generator(None).unwrap();
        let base = ERROR_THRESHOLDS[usize::from(b'+')]; // Q10
        let qualities = vec![b'+'; 100]; // 10 expected errors

        let bases = generator.error_model.mean_error_bases();
        // Without runs, errors are drawn so that their bases (indels count by length) add up to 10
        let t = generator.scaled_thresholds(&qualities, &[], 10.0, 0.0)[usize::from(b'+')];
        assert!((f64::from(t) / f64::from(base) - 1.0 / bases).abs() < 1e-6);
        // 4 bases expected from run errors leave 6 for the per-base model
        let t = generator.scaled_thresholds(&qualities, &[], 10.0, 4.0)[usize::from(b'+')];
        assert!((f64::from(t) / f64::from(base) - 0.6 / bases).abs() < 1e-6);
        // Inside runs only substitutions apply, so the other bases make up for it
        let t = generator.scaled_thresholds(&qualities, &[0..50], 10.0, 0.0)[usize::from(b'+')];
        let model = &generator.error_model;
        let expected = 1.0 / (0.5 * model.mean_error_bases() + 0.5 * model.substitution_rate);
        assert!((f64::from(t) / f64::from(base) - expected).abs() < 1e-6);
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

        let (_, _, runs, _) =
            generator.apply_homopolymer_errors(seq.clone(), vec![b'I'; seq.len()], 0.0);
        assert_eq!(runs, vec![0..3, 62..66, 84..87]);

        // At Q0 only those runs change, and each changes at least once in a few tries
        let mut changed = std::collections::BTreeSet::new();
        for _ in 0..20 {
            let (out, _, runs, _) = generator.apply_homopolymer_errors(
                seq.clone(),
                vec![b'!'; seq.len()],
                seq.len() as f64,
            );
            let (before, after) = (rle(&seq), rle(&out));
            assert_eq!(before.len(), after.len());
            for (b, a) in before.iter().zip(&after).filter(|(b, a)| b != a) {
                changed.insert(*b);
                assert!(
                    runs.iter()
                        .any(|r| r.len() == a.1 && out[r.clone()].iter().all(|&x| x == a.0))
                );
            }
        }
        assert_eq!(changed, [(b'C', 3), (b'G', 3), (b'T', 4)].into());
    }

    #[test]
    fn test_reference_normalized_to_uppercase_acgt() {
        let mut generator = error_free_generator(b"acgtNNRYacgt", 12, false, 0.0);
        let read = generator.generate_read().unwrap();
        assert!(read.sequence().iter().all(|b| b"ACGT".contains(b)));
    }
}
