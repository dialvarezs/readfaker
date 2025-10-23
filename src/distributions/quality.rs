use rand::Rng;
use std::collections::HashMap;

/// Empirical distribution of quality scores grouped by read length.
pub struct QualityDistribution {
    qualities_by_length: HashMap<usize, Vec<Vec<u8>>>,
}

impl Default for QualityDistribution {
    fn default() -> Self {
        Self {
            qualities_by_length: HashMap::new(),
        }
    }
}

impl QualityDistribution {
    /// Creates a new empty quality distribution.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds observed quality scores to the distribution.
    ///
    /// # Arguments
    /// * `length` - Length of the read
    /// * `quality` - Quality score string for this read
    pub fn add_value(&mut self, length: usize, quality: Vec<u8>) {
        self.qualities_by_length
            .entry(length)
            .or_insert(Vec::new())
            .push(quality);
    }

    /// Samples quality scores for a given read length.
    ///
    /// Tries to find an exact length match first, then falls back to the
    /// closest available length if no exact match exists.
    ///
    /// # Arguments
    /// * `length` - Desired read length
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// Quality score string, or None if distribution is empty
    pub fn sample<R: Rng>(&self, length: usize, rng: &mut R) -> Option<Vec<u8>> {
        // Try exact match first
        if let Some(qualities) = self.qualities_by_length.get(&length) {
            let index = rng.random_range(0..qualities.len());
            return Some(qualities[index].clone());
        }

        // Find closest length if no exact match
        let closest_length = self
            .qualities_by_length
            .keys()
            .min_by_key(|&&len| (len as i64 - length as i64).abs())?;

        let qualities = &self.qualities_by_length[closest_length];
        let index = rng.random_range(0..qualities.len());
        Some(qualities[index].clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_add_and_sample() {
        let mut dist = QualityDistribution::new();
        dist.add_value(5, vec![b'?'; 5]);

        let mut rng = StdRng::seed_from_u64(42);
        let sampled = dist.sample(5, &mut rng);

        assert!(sampled.is_some());
        let sampled = sampled.unwrap();
        assert_eq!(sampled.len(), 5);
        assert!(sampled.iter().all(|&q| q >= 33));
    }
}
