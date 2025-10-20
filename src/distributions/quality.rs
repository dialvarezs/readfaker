use rand::Rng;
use std::collections::HashMap;

/// Distribution of read qualities sampled from observed data.
pub struct QualityDistribution {
    qualities_by_length: HashMap<usize, Vec<Vec<u8>>>,
}

impl QualityDistribution {
    /// Creates a new empty quality distribution.
    pub fn new() -> Self {
        Self {
            qualities_by_length: HashMap::new(),
        }
    }

    /// Adds observed qualities to the distribution.
    pub fn add_value(&mut self, length: usize, quality: Vec<u8>) {
        self.qualities_by_length
            .entry(length)
            .or_insert(Vec::new())
            .push(quality);
    }

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
