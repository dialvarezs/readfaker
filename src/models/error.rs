use anyhow::{Result, bail};
use rand::{Rng, RngExt};

// Defaults measured on R10.4.1 E8.2 (dorado fast/hac/sup) reads aligned to their isolate assemblies;
// the error-event mix was stable across basecaller tiers.
const SUBSTITUTION_DEFAULT_RATE: f64 = 0.55;
const INSERTION_DEFAULT_RATE: f64 = 0.18;
const DELETION_DEFAULT_RATE: f64 = 0.27;
const INSERTION_EXTENSION_DEFAULT_RATE: f64 = 0.27;
const DELETION_EXTENSION_DEFAULT_RATE: f64 = 0.3;
const EXTENSION_LIMIT: usize = 100;

/// Extra homopolymer length-error weight by run length (index = length, last entry for longer runs).
/// A run of length L gets an extra length error with probability `1 - (1 - p)^w(L)`, where `p` is the
/// read's mean per-base error probability. Weights are fitted so that, added to the per-base errors,
/// homopolymer accuracy matches real R10.4.1 sup/hac reads from bacterial isolates (runs up to 9)
/// and HG002 against its Q100 assembly (runs up to 12+); geometric mean across datasets and tiers.
const HOMOPOLYMER_WEIGHTS: [f64; 13] = [
    0.0, 0.0, 0.0, 0.4, 1.2, 2.2, 3.7, 5.2, 13.0, 35.0, 50.0, 100.0, 270.0,
];
/// Fraction of homopolymer length errors that shorten the run (the rest lengthen it).
const HOMOPOLYMER_SHORTEN_RATE: f64 = 0.75;
/// Long runs are lengthened more often than short ones.
const LONG_HOMOPOLYMER_SHORTEN_RATE: f64 = 0.65;
const LONG_HOMOPOLYMER_LENGTH: usize = 8;

/// Type of sequencing error alteration applied to a nucleotide position.
#[derive(Debug)]
pub enum AlterationType {
    /// Replace nucleotide with a different one
    Substitution,
    /// Insert random nucleotides after the current position
    Insertion(usize),
    /// Delete nucleotides starting at the current position
    Deletion(usize),
}

/// Model defining the relative probabilities of different sequencing error types.
///
/// # Rate Constraints
/// - All rates must be in the range [0.0, 1.0]
/// - The sum of substitution, insertion, and deletion rates must be ≤ 1.0
/// - Extension rates define geometric distribution for indel lengths
/// - If the sum < 1.0, some errors will result in no alteration
#[derive(Debug)]
pub struct ErrorModel {
    pub substitution_rate: f64,
    pub insertion_rate: f64,
    pub deletion_rate: f64,
    pub insertion_extension_rate: f64,
    pub deletion_extension_rate: f64,
}

impl ErrorModel {
    /// Creates a new error model with specified or default rates.
    ///
    /// # Arguments
    /// * `substitution_rate` - Probability of substitution errors (default: 0.55)
    /// * `insertion_rate` - Probability of insertion errors (default: 0.18)
    /// * `deletion_rate` - Probability of deletion errors (default: 0.27)
    /// * `insertion_extension_rate` - Probability of extending an insertion (default: 0.27)
    /// * `deletion_extension_rate` - Probability of extending a deletion (default: 0.3)
    ///
    /// # Returns
    /// A validated `ErrorModel` instance
    ///
    /// # Errors
    /// Returns an error if:
    /// - Any rate is outside [0.0, 1.0]
    /// - The sum of error rates exceeds 1.0
    ///
    /// # Example
    /// ```
    /// use readfaker::models::ErrorModel;
    ///
    /// // Use default rates (0.55, 0.18, 0.27, 0.27, 0.3)
    /// let model = ErrorModel::new(None, None, None, None, None).unwrap();
    ///
    /// // Custom rates with longer indels
    /// let model = ErrorModel::new(
    ///     Some(0.5), Some(0.3), Some(0.1),
    ///     Some(0.25), Some(0.25)
    /// ).unwrap();
    /// ```
    pub fn new(
        substitution_rate: Option<f64>,
        insertion_rate: Option<f64>,
        deletion_rate: Option<f64>,
        insertion_extension_rate: Option<f64>,
        deletion_extension_rate: Option<f64>,
    ) -> Result<Self> {
        let substitution = substitution_rate.unwrap_or(SUBSTITUTION_DEFAULT_RATE);
        let insertion = insertion_rate.unwrap_or(INSERTION_DEFAULT_RATE);
        let deletion = deletion_rate.unwrap_or(DELETION_DEFAULT_RATE);
        let ins_ext = insertion_extension_rate.unwrap_or(INSERTION_EXTENSION_DEFAULT_RATE);
        let del_ext = deletion_extension_rate.unwrap_or(DELETION_EXTENSION_DEFAULT_RATE);

        for (name, rate) in [
            ("Substitution rate", substitution),
            ("Insertion rate", insertion),
            ("Deletion rate", deletion),
            ("Insertion extension rate", ins_ext),
            ("Deletion extension rate", del_ext),
        ] {
            if !(0.0..=1.0).contains(&rate) {
                bail!("{} must be between 0.0 and 1.0, got {}", name, rate);
            }
        }

        let sum = substitution + insertion + deletion;
        if sum > 1.0 {
            bail!(
                "Error rates must sum to at most 1.0 (got substitution={}, insertion={}, deletion={}, sum={})",
                substitution,
                insertion,
                deletion,
                sum
            );
        }

        Ok(ErrorModel {
            substitution_rate: substitution,
            insertion_rate: insertion,
            deletion_rate: deletion,
            insertion_extension_rate: ins_ext,
            deletion_extension_rate: del_ext,
        })
    }

    /// Sample length from geometric distribution based on extension rate
    fn sample_length(&self, rng: &mut impl Rng, extension_rate: f64) -> usize {
        let mut length = 1;
        // Simple geometric sampling: keep extending while random < rate
        // If rate is 0.0, this loop never runs, returns 1
        while rng.random_range(0.0..1.0) < extension_rate && length < EXTENSION_LIMIT {
            length += 1;
        }
        length
    }

    /// Computes the probability of a homopolymer length error for each run length.
    ///
    /// Computed once per read, so runs can be checked with a table lookup.
    ///
    /// # Arguments
    /// * `error_probability` - Mean per-base error probability of the read
    ///
    /// # Returns
    /// Array indexed by run length (the last entry applies to longer runs); 0.0 for runs shorter
    /// than 3, which the per-base model covers
    pub fn homopolymer_error_probabilities(
        &self,
        error_probability: f64,
    ) -> [f64; HOMOPOLYMER_WEIGHTS.len()] {
        HOMOPOLYMER_WEIGHTS.map(|weight| 1.0 - (1.0 - error_probability).powf(weight))
    }

    /// Samples a length change for a homopolymer run that has an error.
    ///
    /// Short runs are mostly shortened, long runs (8+) nearly as often lengthened. The size of the
    /// change follows the same geometric distributions as regular indels.
    ///
    /// # Arguments
    /// * `length` - Length of the homopolymer run
    /// * `rng` - Random number generator for sampling
    ///
    /// # Returns
    /// Signed change in run length: negative shortens the run, positive lengthens it
    pub fn get_homopolymer_alteration(&self, length: usize, rng: &mut impl Rng) -> isize {
        let shorten_rate = if length >= LONG_HOMOPOLYMER_LENGTH {
            LONG_HOMOPOLYMER_SHORTEN_RATE
        } else {
            HOMOPOLYMER_SHORTEN_RATE
        };
        if rng.random_bool(shorten_rate) {
            -(self.sample_length(rng, self.deletion_extension_rate) as isize)
        } else {
            self.sample_length(rng, self.insertion_extension_rate) as isize
        }
    }

    /// Randomly determines which type of error alteration to apply.
    ///
    /// Uses the model's rate probabilities to select between substitution, insertion,
    /// and deletion. If the sum of rates is less than 1.0, this may return `None`,
    /// indicating no alteration should be applied.
    ///
    /// # Arguments
    /// * `rng` - Random number generator for sampling
    ///
    /// # Returns
    /// * `Some(AlterationType)` - The type of alteration to apply
    /// * `None` - No alteration (when random value exceeds the sum of all rates)
    pub fn get_alteration_type(&self, rng: &mut impl Rng) -> Option<AlterationType> {
        let r = rng.random_range(0.0..1.0);

        if r < self.substitution_rate {
            Some(AlterationType::Substitution)
        } else if r < self.substitution_rate + self.insertion_rate {
            let len = self.sample_length(rng, self.insertion_extension_rate);
            Some(AlterationType::Insertion(len))
        } else if r < self.substitution_rate + self.insertion_rate + self.deletion_rate {
            let len = self.sample_length(rng, self.deletion_extension_rate);
            Some(AlterationType::Deletion(len))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_custom_valid_values() {
        let model = ErrorModel::new(Some(0.5), Some(0.3), Some(0.2), Some(0.5), Some(0.5)).unwrap();
        assert_eq!(model.substitution_rate, 0.5);
        assert_eq!(model.insertion_rate, 0.3);
        assert_eq!(model.deletion_rate, 0.2);
        assert_eq!(model.insertion_extension_rate, 0.5);
        assert_eq!(model.deletion_extension_rate, 0.5);
    }

    #[test]
    fn test_sum_exceeds_one() {
        let result = ErrorModel::new(Some(0.6), Some(0.3), Some(0.2), None, None);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("must sum to at most 1.0"));
    }

    #[test]
    fn test_get_alteration_type_substitution() {
        let model = ErrorModel::new(Some(1.0), Some(0.0), Some(0.0), None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        // With 100% substitution rate, should always return Substitution
        for _ in 0..10 {
            let alteration = model.get_alteration_type(&mut rng);
            assert!(matches!(alteration, Some(AlterationType::Substitution)));
        }
    }

    #[test]
    fn test_get_alteration_type_insertion() {
        let model = ErrorModel::new(Some(0.0), Some(1.0), Some(0.0), None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        // With 100% insertion rate, should always return Insertion
        for _ in 0..10 {
            let alteration = model.get_alteration_type(&mut rng);
            match alteration {
                Some(AlterationType::Insertion(count)) => {
                    assert!(count >= 1); // Default extension is 0.4, so length >= 1
                }
                _ => panic!("Expected Insertion, got {:?}", alteration),
            }
        }
    }

    #[test]
    fn test_get_alteration_type_deletion() {
        let model = ErrorModel::new(Some(0.0), Some(0.0), Some(1.0), None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        // With 100% deletion rate, should always return Deletion
        for _ in 0..10 {
            let alteration = model.get_alteration_type(&mut rng);
            match alteration {
                Some(AlterationType::Deletion(count)) => {
                    assert!(count >= 1); // Default extension is 0.4, so length >= 1
                }
                _ => panic!("Expected Deletion, got {:?}", alteration),
            }
        }
    }

    #[test]
    fn test_get_alteration_type_mixed() {
        let model = ErrorModel::new(Some(0.5), Some(0.3), Some(0.2), None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        let mut substitution_count = 0;
        let mut insertion_count = 0;
        let mut deletion_count = 0;
        let mut none_count = 0;

        // Sample many times to check distribution
        for _ in 0..1000 {
            match model.get_alteration_type(&mut rng) {
                Some(AlterationType::Substitution) => substitution_count += 1,
                Some(AlterationType::Insertion(_)) => insertion_count += 1,
                Some(AlterationType::Deletion(_)) => deletion_count += 1,
                None => none_count += 1,
            }
        }

        // Should have substitutions (most common), insertions, and deletions
        assert!(substitution_count > 400); // Expect ~500
        assert!(insertion_count > 200); // Expect ~300
        assert!(deletion_count > 100); // Expect ~200
        assert_eq!(none_count, 0); // Sum = 1.0, so no None
    }

    #[test]
    fn test_indel_extension() {
        // Set high extension rate to ensure we get length > 1
        let model = ErrorModel::new(Some(0.0), Some(0.5), Some(0.5), Some(0.8), Some(0.8)).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        let mut max_ins_len = 0;
        let mut max_del_len = 0;

        for _ in 0..1000 {
            match model.get_alteration_type(&mut rng) {
                Some(AlterationType::Insertion(count)) => max_ins_len = max_ins_len.max(count),
                Some(AlterationType::Deletion(count)) => max_del_len = max_del_len.max(count),
                _ => {}
            }
        }

        assert!(max_ins_len > 1);
        assert!(max_del_len > 1);
    }

    #[test]
    fn test_homopolymer_error_probabilities_grow_with_length() {
        let model = ErrorModel::new(None, None, None, None, None).unwrap();
        let probabilities = model.homopolymer_error_probabilities(0.01);
        assert_eq!(probabilities[..3], [0.0; 3]);
        assert!(probabilities[3..].windows(2).all(|w| w[0] < w[1]));
        assert!(probabilities.iter().all(|&p| (0.0..=1.0).contains(&p)));
        assert_eq!(model.homopolymer_error_probabilities(1.0)[3], 1.0);
    }

    #[test]
    fn test_indel_length_capped_at_limit() {
        let model = ErrorModel::new(None, None, None, None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        assert_eq!(model.sample_length(&mut rng, 1.0), EXTENSION_LIMIT);
    }
}
