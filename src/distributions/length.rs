use rand::Rng;
use std::collections::HashMap;

/// Distribution of read lengths sampled from observed data.
pub struct LengthDistribution {
    length_histogram: HashMap<usize, usize>,
    total_count: usize,
}

impl LengthDistribution {
    /// Creates a new empty length distribution.
    pub fn new() -> Self {
        Self {
            length_histogram: HashMap::new(),
            total_count: 0,
        }
    }

    /// Adds observed lengths to the distribution.
    pub fn add_value(&mut self, length: usize) {
        self.length_histogram
            .entry(length)
            .and_modify(|c| *c += 1)
            .or_insert(1);
        self.total_count += 1;
    }

    /// Samples a random length from the distribution.
    pub fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        let target = rng.random_range(0..self.total_count);
        let mut cumulative = 0;

        for (&length, &count) in &self.length_histogram {
            cumulative += count;
            if cumulative > target {
                return length;
            }
        }

        unreachable!();
    }
}
