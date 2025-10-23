use rand::Rng;
use std::collections::HashMap;

/// Empirical model of quality scores built from observed reads, grouped by read length.
#[derive(Default)]
pub struct QualityModel {
    qualities_by_length: HashMap<usize, Vec<Vec<u8>>>,
}

impl QualityModel {
    /// Creates a new empty quality model.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds observed quality scores to the empirical model.
    ///
    /// # Arguments
    /// * `length` - Length of the read
    /// * `quality` - Quality score string for this read
    pub fn add_value(&mut self, length: usize, quality: Vec<u8>) {
        self.qualities_by_length
            .entry(length)
            .or_default()
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
    /// Quality score string, or None if model is empty
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
            .min_by_key(|&&len| len.abs_diff(length))?;

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
        let mut dist = QualityModel::new();
        dist.add_value(5, vec![b'?'; 5]);

        let mut rng = StdRng::seed_from_u64(42);
        let sampled = dist.sample(5, &mut rng);

        assert!(sampled.is_some());
        let sampled = sampled.unwrap();
        assert_eq!(sampled.len(), 5);
        assert!(sampled.iter().all(|&q| q >= 33));
    }

    #[test]
    fn test_sample_fallback_to_closest() {
        // Test that sampling falls back to closest length when exact match not found
        let mut dist = QualityModel::new();
        dist.add_value(100, vec![b'I'; 100]); // Length 100
        dist.add_value(200, vec![b'J'; 200]); // Length 200
        dist.add_value(500, vec![b'K'; 500]); // Length 500

        let mut rng = StdRng::seed_from_u64(42);

        // Request length 150 - closest is 100 or 200 (both 50 away, min_by_key picks first)
        let sampled = dist.sample(150, &mut rng);
        assert!(sampled.is_some());
        let sampled = sampled.unwrap();
        // Should get quality string from either 100 or 200
        assert!(sampled.len() == 100 || sampled.len() == 200);

        // Request length 450 - closest is 500 (50 away)
        let sampled = dist.sample(450, &mut rng);
        assert!(sampled.is_some());
        assert_eq!(sampled.unwrap().len(), 500);
    }
}
