use anyhow::{bail, Result};
use rand::Rng;

const SUBSTITUTION_DEFAULT_RATE: f64 = 0.7;
const INSERTION_DEFAULT_RATE: f64 = 0.1;
const DELETION_DEFAULT_RATE: f64 = 0.2;

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
/// - The sum of all rates must be ≤ 1.0
/// - If the sum < 1.0, some errors will result in no alteration
#[derive(Debug)]
pub struct ErrorModel {
    pub substitution_rate: f64,
    pub insertion_rate: f64,
    pub deletion_rate: f64,
}

impl ErrorModel {
    /// Creates a new error model with specified or default rates.
    ///
    /// # Arguments
    /// * `substitution_rate` - Probability of substitution errors (default: 0.7)
    /// * `insertion_rate` - Probability of insertion errors (default: 0.1)
    /// * `deletion_rate` - Probability of deletion errors (default: 0.2)
    ///
    /// # Returns
    /// A validated `ErrorModel` instance
    ///
    /// # Errors
    /// Returns an error if:
    /// - Any rate is outside [0.0, 1.0]
    /// - The sum of all rates exceeds 1.0
    ///
    /// # Example
    /// ```
    /// use readfaker::models::ErrorModel;
    ///
    /// // Use default rates (0.7, 0.1, 0.2)
    /// let model = ErrorModel::new(None, None, None).unwrap();
    ///
    /// // Custom rates
    /// let model = ErrorModel::new(Some(0.5), Some(0.3), Some(0.1)).unwrap();
    /// ```
    pub fn new(
        substitution_rate: Option<f64>,
        insertion_rate: Option<f64>,
        deletion_rate: Option<f64>,
    ) -> Result<Self> {
        let substitution = substitution_rate.unwrap_or(SUBSTITUTION_DEFAULT_RATE);
        let insertion = insertion_rate.unwrap_or(INSERTION_DEFAULT_RATE);
        let deletion = deletion_rate.unwrap_or(DELETION_DEFAULT_RATE);

        if !(0.0..=1.0).contains(&substitution) {
            bail!(
                "Substitution rate must be between 0.0 and 1.0, got {}",
                substitution
            );
        }
        if !(0.0..=1.0).contains(&insertion) {
            bail!(
                "Insertion rate must be between 0.0 and 1.0, got {}",
                insertion
            );
        }
        if !(0.0..=1.0).contains(&deletion) {
            bail!(
                "Deletion rate must be between 0.0 and 1.0, got {}",
                deletion
            );
        }

        let sum = substitution + insertion + deletion;
        if sum > 1.0 {
            bail!(
                "Error rates must sum to at most 1.0 (got substitution={}, insertion={}, deletion={}, sum={})",
                substitution, insertion, deletion, sum
            );
        }

        Ok(ErrorModel {
            substitution_rate: substitution,
            insertion_rate: insertion,
            deletion_rate: deletion,
        })
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
    ///
    /// # Current Implementation
    /// Always returns count = 1 for insertions and deletions. Multi-nucleotide
    /// operations will be supported in future versions.
    pub fn get_alteration_type(&self, rng: &mut impl Rng) -> Option<AlterationType> {
        let r = rng.random_range(0.0..1.0);

        if r < self.substitution_rate {
            Some(AlterationType::Substitution)
        } else if r < self.substitution_rate + self.insertion_rate {
            Some(AlterationType::Insertion(1))
        } else if r < self.substitution_rate + self.insertion_rate + self.deletion_rate {
            Some(AlterationType::Deletion(1))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_default_values() {
        let model = ErrorModel::new(None, None, None).unwrap();
        assert_eq!(model.substitution_rate, 0.7);
        assert_eq!(model.insertion_rate, 0.1);
        assert_eq!(model.deletion_rate, 0.2);
    }

    #[test]
    fn test_custom_valid_values() {
        let model = ErrorModel::new(Some(0.5), Some(0.3), Some(0.2)).unwrap();
        assert_eq!(model.substitution_rate, 0.5);
        assert_eq!(model.insertion_rate, 0.3);
        assert_eq!(model.deletion_rate, 0.2);
    }

    #[test]
    fn test_sum_exceeds_one() {
        let result = ErrorModel::new(Some(0.6), Some(0.3), Some(0.2));
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("must sum to at most 1.0"));
        assert!(err_msg.contains("substitution=0.6"));
        assert!(err_msg.contains("insertion=0.3"));
        assert!(err_msg.contains("deletion=0.2"));
    }

    #[test]
    fn test_get_alteration_type_substitution() {
        let model = ErrorModel::new(Some(1.0), Some(0.0), Some(0.0)).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        // With 100% substitution rate, should always return Substitution
        for _ in 0..10 {
            let alteration = model.get_alteration_type(&mut rng);
            assert!(matches!(alteration, Some(AlterationType::Substitution)));
        }
    }

    #[test]
    fn test_get_alteration_type_insertion() {
        let model = ErrorModel::new(Some(0.0), Some(1.0), Some(0.0)).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        // With 100% insertion rate, should always return Insertion
        for _ in 0..10 {
            let alteration = model.get_alteration_type(&mut rng);
            match alteration {
                Some(AlterationType::Insertion(count)) => {
                    assert_eq!(count, 1);
                }
                _ => panic!("Expected Insertion, got {:?}", alteration),
            }
        }
    }

    #[test]
    fn test_get_alteration_type_deletion() {
        let model = ErrorModel::new(Some(0.0), Some(0.0), Some(1.0)).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        // With 100% deletion rate, should always return Deletion
        for _ in 0..10 {
            let alteration = model.get_alteration_type(&mut rng);
            match alteration {
                Some(AlterationType::Deletion(count)) => {
                    assert_eq!(count, 1);
                }
                _ => panic!("Expected Deletion, got {:?}", alteration),
            }
        }
    }

    #[test]
    fn test_get_alteration_type_mixed() {
        let model = ErrorModel::new(Some(0.5), Some(0.3), Some(0.2)).unwrap();
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
}
