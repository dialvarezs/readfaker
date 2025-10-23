use rand::Rng;
use std::collections::HashMap;

/// Empirical distribution of read lengths built from observed data.
pub struct LengthDistribution {
    length_histogram: HashMap<usize, usize>,
    total_count: usize,
}

impl Default for LengthDistribution {
    fn default() -> Self {
        Self {
            length_histogram: HashMap::new(),
            total_count: 0,
        }
    }
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
}
