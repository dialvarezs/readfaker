use anyhow::{Result, bail};
use rand::{Rng, RngExt};

// Defaults measured on R10.4.1 E8.2 reads (dorado fast/hac/sup) from bacterial isolates aligned to
// their assemblies and HG002 aligned to its Q100 assembly, counting only errors outside
// homopolymer runs (those follow the homopolymer model); mean across datasets.
const SUBSTITUTION_DEFAULT_RATE: f64 = 0.59;
const INSERTION_DEFAULT_RATE: f64 = 0.18;
const DELETION_DEFAULT_RATE: f64 = 0.23;
const INSERTION_EXTENSION_DEFAULT_RATE: f64 = 0.26;
const DELETION_EXTENSION_DEFAULT_RATE: f64 = 0.34;
const EXTENSION_LIMIT: usize = 100;

/// Homopolymer length-error model: `(a, b)` by run length (index = length, last entry for longer
/// runs). A run has a length error with probability `1 - exp(-a * p^b)`, where `p` is the read's
/// mean per-base error probability. `b` falls with run length because errors in long runs depend
/// little on read quality. Runs shorter than 3 are left to the per-base model.
///
/// Fitted on R10.4.1 reads (dorado fast/hac/sup) from bacterial isolates and HG002.
const HOMOPOLYMER_MODEL: [(f64, f64); 13] = [
    (0.0, 1.0),
    (0.0, 1.0),
    (0.0, 1.0),
    (1.12, 0.93),
    (1.58, 0.87),
    (1.76, 0.80),
    (2.06, 0.73),
    (2.04, 0.66),
    (2.13, 0.53),
    (2.48, 0.42),
    (2.28, 0.32),
    (1.94, 0.19),
    (1.89, 0.058),
];
/// Fraction of homopolymer length errors that shorten the run (the rest lengthen it), by run length
/// (index = length, last entry for longer runs).
const HOMOPOLYMER_SHORTEN_RATES: [f64; 7] = [0.0, 0.0, 0.0, 0.50, 0.60, 0.65, 0.70];
/// Probability of shortening or lengthening a run by one more base.
const HOMOPOLYMER_SHORTEN_EXTENSION_RATE: f64 = 0.2;
const HOMOPOLYMER_LENGTHEN_EXTENSION_RATE: f64 = 0.25;

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
/// - Substitution, insertion, and deletion rates must sum to 1.0
/// - Extension rates define geometric distribution for indel lengths
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
    /// * `substitution_rate` - Probability of substitution errors (default: 0.59)
    /// * `insertion_rate` - Probability of insertion errors (default: 0.18)
    /// * `deletion_rate` - Probability of deletion errors (default: 0.23)
    /// * `insertion_extension_rate` - Probability of extending an insertion (default: 0.26)
    /// * `deletion_extension_rate` - Probability of extending a deletion (default: 0.34)
    ///
    /// # Returns
    /// A validated `ErrorModel` instance
    ///
    /// # Errors
    /// Returns an error if:
    /// - Any rate is outside [0.0, 1.0]
    /// - The error rates don't sum to 1.0
    ///
    /// # Example
    /// ```
    /// use readfaker::models::ErrorModel;
    ///
    /// // Use default rates (0.59, 0.18, 0.23, 0.26, 0.34)
    /// let model = ErrorModel::new(None, None, None, None, None).unwrap();
    ///
    /// // Custom rates with longer indels
    /// let model = ErrorModel::new(
    ///     Some(0.5), Some(0.3), Some(0.2),
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
        if (sum - 1.0).abs() > 1e-6 {
            bail!(
                "Error rates must sum to 1.0 (got substitution={}, insertion={}, deletion={}, sum={})",
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

    /// Mean length of indels sampled with `sample_length`.
    fn mean_length(extension_rate: f64) -> f64 {
        if extension_rate >= 1.0 {
            return EXTENSION_LIMIT as f64;
        }
        (1.0 - extension_rate.powi(EXTENSION_LIMIT as i32)) / (1.0 - extension_rate)
    }

    /// Mean number of bases changed by an error from `get_alteration_type`: one for a
    /// substitution, the indel length for insertions and deletions.
    pub fn mean_error_bases(&self) -> f64 {
        self.substitution_rate
            + self.insertion_rate * Self::mean_length(self.insertion_extension_rate)
            + self.deletion_rate * Self::mean_length(self.deletion_extension_rate)
    }

    /// Mean number of bases added or removed by a homopolymer length error.
    ///
    /// # Arguments
    /// * `length` - Length of the homopolymer run
    pub fn mean_homopolymer_error_bases(&self, length: usize) -> f64 {
        let shorten_rate = Self::shorten_rate(length);
        shorten_rate
            * Self::mean_length(HOMOPOLYMER_SHORTEN_EXTENSION_RATE).min((length - 1) as f64)
            + (1.0 - shorten_rate) * Self::mean_length(HOMOPOLYMER_LENGTHEN_EXTENSION_RATE)
    }

    fn shorten_rate(length: usize) -> f64 {
        HOMOPOLYMER_SHORTEN_RATES[length.min(HOMOPOLYMER_SHORTEN_RATES.len() - 1)]
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
    ) -> [f64; HOMOPOLYMER_MODEL.len()] {
        HOMOPOLYMER_MODEL.map(|(a, b)| 1.0 - (-a * error_probability.powf(b)).exp())
    }

    /// Samples a length change for a homopolymer run that has an error.
    ///
    /// Runs of 3 are as often shortened as lengthened; longer runs are mostly shortened. The size of
    /// the change follows a geometric distribution, with lower extension rates than regular indels.
    ///
    /// # Arguments
    /// * `length` - Length of the homopolymer run
    /// * `rng` - Random number generator for sampling
    ///
    /// # Returns
    /// Signed change in run length: negative shortens the run, positive lengthens it
    pub fn get_homopolymer_alteration(&self, length: usize, rng: &mut impl Rng) -> isize {
        if rng.random_bool(Self::shorten_rate(length)) {
            -(self.sample_length(rng, HOMOPOLYMER_SHORTEN_EXTENSION_RATE) as isize)
        } else {
            self.sample_length(rng, HOMOPOLYMER_LENGTHEN_EXTENSION_RATE) as isize
        }
    }

    /// Randomly determines which type of error alteration to apply.
    ///
    /// Uses the model's rate probabilities to select between substitution, insertion,
    /// and deletion.
    ///
    /// # Arguments
    /// * `rng` - Random number generator for sampling
    ///
    /// # Returns
    /// The type of alteration to apply
    pub fn get_alteration_type(&self, rng: &mut impl Rng) -> AlterationType {
        let r = rng.random_range(0.0..1.0);

        if r < self.substitution_rate {
            AlterationType::Substitution
        } else if r < self.substitution_rate + self.insertion_rate {
            let len = self.sample_length(rng, self.insertion_extension_rate);
            AlterationType::Insertion(len)
        } else {
            let len = self.sample_length(rng, self.deletion_extension_rate);
            AlterationType::Deletion(len)
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
    fn test_sum_must_be_one() {
        for (sub, ins, del) in [(0.6, 0.3, 0.2), (0.5, 0.2, 0.2)] {
            let result = ErrorModel::new(Some(sub), Some(ins), Some(del), None, None);
            let err_msg = result.unwrap_err().to_string();
            assert!(err_msg.contains("must sum to 1.0"));
        }
    }

    #[test]
    fn test_get_alteration_type_substitution() {
        let model = ErrorModel::new(Some(1.0), Some(0.0), Some(0.0), None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        // With 100% substitution rate, should always return Substitution
        for _ in 0..10 {
            let alteration = model.get_alteration_type(&mut rng);
            assert!(matches!(alteration, AlterationType::Substitution));
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
                AlterationType::Insertion(count) => {
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
                AlterationType::Deletion(count) => {
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

        // Sample many times to check distribution
        for _ in 0..1000 {
            match model.get_alteration_type(&mut rng) {
                AlterationType::Substitution => substitution_count += 1,
                AlterationType::Insertion(_) => insertion_count += 1,
                AlterationType::Deletion(_) => deletion_count += 1,
            }
        }

        // Should have substitutions (most common), insertions, and deletions
        assert!(substitution_count > 400); // Expect ~500
        assert!(insertion_count > 200); // Expect ~300
        assert!(deletion_count > 100); // Expect ~200
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
                AlterationType::Insertion(count) => max_ins_len = max_ins_len.max(count),
                AlterationType::Deletion(count) => max_del_len = max_del_len.max(count),
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
        // Long runs depend less on read quality than short ones
        let low = model.homopolymer_error_probabilities(0.1);
        assert!(low[3] / probabilities[3] > low[9] / probabilities[9]);
    }

    #[test]
    fn test_mean_homopolymer_error_bases_matches_sampling() {
        let model = ErrorModel::new(None, None, None, None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(7);
        let n = 200_000;
        let changed: usize = (0..n)
            .map(|_| {
                model
                    .get_homopolymer_alteration(6, &mut rng)
                    .unsigned_abs()
                    .min(5)
            })
            .sum();
        let mean = changed as f64 / n as f64;
        // capping a shortening at 5 bases barely matters at these extension rates
        assert!(
            (mean - model.mean_homopolymer_error_bases(6)).abs() < 0.01,
            "{mean}"
        );
    }

    #[test]
    fn test_indel_length_capped_at_limit() {
        let model = ErrorModel::new(None, None, None, None, None).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        assert_eq!(model.sample_length(&mut rng, 1.0), EXTENSION_LIMIT);
    }
}
