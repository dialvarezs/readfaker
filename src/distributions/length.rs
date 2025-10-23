use rand::Rng;
use std::collections::BTreeMap;

/// Empirical distribution of read lengths built from observed data.
#[derive(Default)]
pub struct LengthDistribution {
    length_histogram: BTreeMap<usize, usize>,
    total_count: usize,
}

impl LengthDistribution {
    /// Creates a new empty length distribution.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an observed length to the distribution.
    ///
    /// # Arguments
    /// * `length` - Read length to add
    pub fn add_value(&mut self, length: usize) {
        self.length_histogram
            .entry(length)
            .and_modify(|c| *c += 1)
            .or_insert(1);
        self.total_count += 1;
    }

    /// Samples a random length from the distribution.
    ///
    /// # Arguments
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// A randomly sampled read length, or None if the distribution is empty
    pub fn sample<R: Rng>(&self, rng: &mut R) -> Option<usize> {
        if self.total_count == 0 {
            return None;
        }

        let target = rng.random_range(0..self.total_count);
        let mut cumulative = 0;

        for (&length, &count) in &self.length_histogram {
            cumulative += count;
            if cumulative > target {
                return Some(length);
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_add_and_sample() {
        let mut dist = LengthDistribution::new();
        dist.add_value(100);
        dist.add_value(100);
        dist.add_value(200);

        let mut rng = StdRng::seed_from_u64(42);
        let sampled = dist.sample(&mut rng).unwrap();
        assert!(sampled == 100 || sampled == 200);
    }

    #[test]
    fn test_deterministic_sampling() {
        // Verify that sampling is reproducible with same seed
        let mut dist = LengthDistribution::new();
        dist.add_value(50);
        dist.add_value(100);
        dist.add_value(150);
        dist.add_value(200);
        dist.add_value(250);

        // Sample with first RNG
        let mut rng1 = StdRng::seed_from_u64(12345);
        let samples1: Vec<usize> = (0..10).map(|_| dist.sample(&mut rng1).unwrap()).collect();

        // Sample with second RNG (same seed)
        let mut rng2 = StdRng::seed_from_u64(12345);
        let samples2: Vec<usize> = (0..10).map(|_| dist.sample(&mut rng2).unwrap()).collect();

        // Should produce identical sequences
        assert_eq!(samples1, samples2);
    }
}
