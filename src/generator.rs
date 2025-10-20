use crate::distributions::{LengthDistribution, QualityDistribution};
use crate::io::fasta::FastaRecord;
use crate::io::fastq::FastqReader;
use crate::io::FastqRecord;
use crate::utils::{get_random_nucleotide, QUALITY_MAPPING};
use anyhow::{anyhow, Result};
use rand::prelude::IndexedRandom;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::path::PathBuf;

/// Generates synthetic reads from reference sequences with realistic quality profiles.
///
/// Randomly samples subsequences from the provided reference sequences, applies
/// sequencing errors based on quality scores, and generates FASTQ records.
///
/// # Arguments
/// * `sequences` - Reference sequences to sample from
/// * `length_distribution` - Distribution for sampling read lengths
/// * `quality_distribution` - Distribution for sampling quality scores
/// * `number_of_reads` - Number of reads to generate
/// * `seed` - Optional random seed for reproducibility (defaults to 0)
///
/// # Returns
/// Vector of generated FASTQ records with simulated errors
pub fn generate_reads(
    sequences: Vec<FastaRecord>,
    length_distribution: LengthDistribution,
    quality_distribution: QualityDistribution,
    number_of_reads: usize,
    seed: Option<u64>,
) -> Result<Vec<FastqRecord>> {
    let mut rng = StdRng::seed_from_u64(seed.unwrap_or(0));
    let mut records = Vec::new();
    let mut read_count = 0;

    while read_count < number_of_reads {
        let length = length_distribution
            .sample(&mut rng)
            .ok_or_else(|| anyhow!("Length distribution is empty"))?;
        let reference_sequence = sequences.choose(&mut rng).unwrap();

        // Skip if sampled length is longer than reference sequence
        if length > reference_sequence.sequence.len() {
            continue;
        }

        let max_start = reference_sequence.sequence.len() - length;
        let start_position = rng.random_range(0..=max_start);
        let mut sequence =
            reference_sequence.sequence[start_position..start_position + length].to_vec();

        let Some(qualities) = quality_distribution.sample(length, &mut rng) else {
            continue; // Skip if no quality string available
        };

        // Insert sequence nucleotide substitutions depending on error value
        qualities.iter().enumerate().for_each(|(i, &quality_ascii)| {
            let phred = quality_ascii.saturating_sub(33);
            if rng.random_range(0.0..1.0) <= QUALITY_MAPPING[&phred] {
                sequence[i] = get_random_nucleotide(sequence[i], &mut rng);
            }
        });

        records.push(FastqRecord {
            id: format!("read_{}", read_count),
            sequence,
            quality: qualities,
        });
        read_count += 1;
    }

    Ok(records)
}

/// Loads length and quality distributions from an existing FASTQ file.
///
/// Reads all records from the input file and builds empirical distributions
/// for read lengths and quality scores.
///
/// # Arguments
/// * `fastq_path` - Path to the FASTQ file to analyze
///
/// # Returns
/// Tuple of (LengthDistribution, QualityDistribution) built from the input file
pub fn load_distributions(
    fastq_path: PathBuf,
) -> Result<(LengthDistribution, QualityDistribution)> {
    let mut length_distribution = LengthDistribution::new();
    let mut quality_distribution = QualityDistribution::new();

    let reader = FastqReader::from_path(&fastq_path)?;

    for record in reader {
        let record = record?;
        length_distribution.add_value(record.len());
        quality_distribution.add_value(record.len(), record.quality);
    }

    Ok((length_distribution, quality_distribution))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_reads() {
        let sequences = vec![FastaRecord {
            id: "seq1".to_string(),
            sequence: b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        }];

        let mut length_dist = LengthDistribution::new();
        let mut quality_dist = QualityDistribution::new();
        length_dist.add_value(10);
        quality_dist.add_value(10, vec![b'?'; 10]); // Phred 30 as ASCII

        let reads = generate_reads(sequences, length_dist, quality_dist, 5, Some(42)).unwrap();

        assert_eq!(reads.len(), 5);
        assert_eq!(reads[0].id, "read_0");
        assert_eq!(reads[0].len(), 10);
        assert!(reads[0].quality.iter().all(|&q| q >= 33));
    }
}
