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
}
