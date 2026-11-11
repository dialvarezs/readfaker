use rand::{Rng, RngExt};
use std::collections::BTreeMap;
use std::sync::OnceLock;

/// Empirical model of read lengths built from observed sequencing data.
#[derive(Default)]
pub struct LengthModel {
    length_histogram: BTreeMap<usize, usize>,
    total_count: usize,
    /// (cumulative count, length) pairs built on first sample, for O(log n) lookups.
    cumulative: OnceLock<Vec<(usize, usize)>>,
}

impl LengthModel {
    /// Creates a new empty length model.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an observed read length to the empirical model.
    ///
    /// # Arguments
    /// * `length` - Read length to add
    pub fn add_value(&mut self, length: usize) {
        *self.length_histogram.entry(length).or_default() += 1;
        self.total_count += 1;
        self.cumulative.take();
    }

    /// Samples a random length from the empirical model.
    ///
    /// # Arguments
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// A randomly sampled read length, or None if the model is empty
    pub fn sample<R: Rng>(&self, rng: &mut R) -> Option<usize> {
        if self.total_count == 0 {
            return None;
        }

        let cumulative = self.cumulative.get_or_init(|| {
            let mut total = 0;
            self.length_histogram
                .iter()
                .map(|(&length, &count)| {
                    total += count;
                    (total, length)
                })
                .collect()
        });

        let target = rng.random_range(0..self.total_count);
        let idx = cumulative.partition_point(|&(total, _)| total <= target);
        Some(cumulative[idx].1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_add_and_sample() {
        let mut model = LengthModel::new();
        model.add_value(100);
        model.add_value(100);
        model.add_value(200);

        let mut rng = StdRng::seed_from_u64(42);
        let sampled = model.sample(&mut rng).unwrap();
        assert!(sampled == 100 || sampled == 200);
    }

    #[test]
    fn test_deterministic_sampling() {
        // Verify that sampling is reproducible with same seed
        let mut model = LengthModel::new();
        model.add_value(50);
        model.add_value(100);
        model.add_value(150);
        model.add_value(200);
        model.add_value(250);

        // Sample with first RNG
        let mut rng1 = StdRng::seed_from_u64(12345);
        let samples1: Vec<usize> = (0..10).map(|_| model.sample(&mut rng1).unwrap()).collect();

        // Sample with second RNG (same seed)
        let mut rng2 = StdRng::seed_from_u64(12345);
        let samples2: Vec<usize> = (0..10).map(|_| model.sample(&mut rng2).unwrap()).collect();

        // Should produce identical sequences
        assert_eq!(samples1, samples2);
    }

    #[test]
    fn test_add_after_sample_updates_distribution() {
        let mut model = LengthModel::new();
        model.add_value(100);
        let mut rng = StdRng::seed_from_u64(42);
        assert_eq!(model.sample(&mut rng), Some(100));

        for _ in 0..1000 {
            model.add_value(200);
        }
        assert!((0..100).any(|_| model.sample(&mut rng) == Some(200)));
    }
}
