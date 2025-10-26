use rand::Rng;
use rand::seq::IndexedRandom;

/// Default maximum read length before grouping into a catch-all bucket (20kb).
const DEFAULT_MAX_BUCKET_LENGTH: usize = 20000;

/// Default width of each length bucket in base pairs (100bp).
const DEFAULT_BUCKET_WIDTH: usize = 100;

/// Default maximum number of quality strings stored per bucket (1000 reads).
const DEFAULT_MAX_ITEMS_PER_BUCKET: usize = 1000;

/// A bucket storing quality strings for reads within a specific length range.
struct QualityBatch {
    qualities: Vec<Vec<u8>>,
    max_capacity: usize,
    total_seen: usize,
}

impl QualityBatch {
    /// Creates a new empty quality batch.
    ///
    /// # Arguments
    /// * `max_capacity` - Maximum number of quality strings to store
    pub fn new(max_capacity: usize) -> Self {
        Self {
            qualities: Vec::new(),
            max_capacity,
            total_seen: 0,
        }
    }

    /// Returns true if the batch has room for more quality strings.
    pub fn has_capacity(&self) -> bool {
        self.qualities.len() < self.max_capacity
    }

    /// Returns an iterator over quality strings that are at least the given length.
    ///
    /// # Arguments
    /// * `length` - Minimum length required
    pub fn qualities_longer_than(&self, length: usize) -> impl Iterator<Item = &Vec<u8>> {
        self.qualities.iter().filter(move |q| q.len() >= length)
    }

    /// Adds a quality string to the batch using reservoir sampling.
    ///
    /// If the batch is not full, the quality string is added directly.
    /// If full, reservoir sampling is used to randomly decide whether to
    /// replace an existing quality string.
    ///
    /// # Arguments
    /// * `quality` - Quality score string as byte vector
    /// * `rng` - Random number generator for reservoir sampling
    pub fn add_value<R: Rng>(&mut self, quality: Vec<u8>, rng: &mut R) {
        self.total_seen += 1;

        if self.has_capacity() {
            self.qualities.push(quality);
        } else {
            // Reservoir sampling: decide whether to include this new quality string
            let random_idx = rng.random_range(0..self.total_seen);
            if random_idx < self.qualities.len() {
                self.qualities[random_idx] = quality;
            }
        }
    }
}

/// Empirical model of quality scores built from observed reads, grouped by read length ranges.
///
/// Quality strings are organized into buckets based on read length to balance memory
/// usage and diversity. Reads are grouped into fixed-width length buckets (default 100bp)
/// up to a threshold (default 20kb), with ultra-long reads stored in a catch-all bucket.
///
/// Each bucket uses reservoir sampling to cap memory usage while maintaining diversity.
pub struct QualityModel {
    /// Quality batches organized by length range.
    batches: Vec<QualityBatch>,
    /// Width of each length bucket in base pairs.
    bucket_width: usize,
    /// Read length threshold for the catch-all bucket.
    max_bucket_length: usize,
}

impl QualityModel {
    /// Creates a new empty quality model with length-based bucketing.
    ///
    /// # Arguments
    /// * `bucket_width` - Width of each length bucket in base pairs (default: 100bp)
    /// * `max_bucket_length` - Maximum length before catch-all bucket (default: 20kb)
    /// * `max_items_per_bucket` - Maximum quality strings per bucket (default: 1000)
    ///
    /// # Example
    /// ```
    /// use readfaker::models::QualityModel;
    ///
    /// // Use defaults (100bp buckets, 20kb threshold, 1000 items per bucket)
    /// let model = QualityModel::new(None, None, None);
    ///
    /// // Custom settings: 200bp buckets, 50kb threshold, 500 items per bucket
    /// let custom_model = QualityModel::new(Some(200), Some(50000), Some(500));
    /// ```
    pub fn new(
        bucket_width: Option<usize>,
        max_bucket_length: Option<usize>,
        max_items_per_bucket: Option<usize>,
    ) -> Self {
        let bucket_width = bucket_width.unwrap_or(DEFAULT_BUCKET_WIDTH);
        let max_bucket_length = max_bucket_length.unwrap_or(DEFAULT_MAX_BUCKET_LENGTH);
        let max_items_per_bucket = max_items_per_bucket.unwrap_or(DEFAULT_MAX_ITEMS_PER_BUCKET);

        let mut batches: Vec<QualityBatch> = (0..max_bucket_length)
            .step_by(bucket_width)
            .map(|_| QualityBatch::new(max_items_per_bucket))
            .collect();
        batches.push(QualityBatch::new(max_items_per_bucket));

        Self {
            batches,
            bucket_width,
            max_bucket_length,
        }
    }

    /// Adds observed quality scores to the empirical model.
    ///
    /// The quality string is assigned to the appropriate length bucket. If the bucket
    /// is full, reservoir sampling determines whether to include it.
    ///
    /// # Arguments
    /// * `length` - Length of the read
    /// * `quality` - Quality score string as byte vector (Phred+33 ASCII encoding)
    /// * `rng` - Random number generator for reservoir sampling
    pub fn add_value<R: Rng>(&mut self, length: usize, quality: Vec<u8>, rng: &mut R) {
        let batch_idx = if length < self.max_bucket_length {
            length / self.bucket_width
        } else {
            self.batches.len() - 1
        };

        self.batches[batch_idx].add_value(quality, rng)
    }

    /// Samples quality scores for a given read length.
    ///
    /// Samples from the appropriate length bucket, filtering for quality strings
    /// that are at least as long as the requested length. If the target bucket is
    /// empty, falls back to the next bucket with longer reads. The sampled quality
    /// string is truncated to exactly match the requested length.
    ///
    /// # Arguments
    /// * `length` - Desired read length
    /// * `rng` - Random number generator
    ///
    /// # Returns
    /// Quality score byte vector truncated to the requested length, or None if no
    /// suitable quality strings are available in any bucket.
    pub fn sample<R: Rng>(&self, length: usize, rng: &mut R) -> Option<Vec<u8>> {
        let mut batch_idx = if length < self.max_bucket_length {
            length / self.bucket_width
        } else {
            self.batches.len() - 1
        };

        loop {
            let candidates: Vec<_> = self.batches[batch_idx]
                .qualities_longer_than(length)
                .collect();

            if !candidates.is_empty() {
                let sampled = candidates.choose(rng).unwrap();
                // Truncate to exact requested length
                return Some(sampled[..length].to_vec());
            }

            // Try next bucket with longer reads
            batch_idx += 1;

            if batch_idx >= self.batches.len() {
                return None;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_bucket_assignment() {
        let mut model = QualityModel::new(None, None, None);
        let mut rng = StdRng::seed_from_u64(42);

        // Add quality strings at different lengths
        model.add_value(50, vec![b'A'; 50], &mut rng);    // Bucket 0 (0-99)
        model.add_value(150, vec![b'B'; 150], &mut rng);  // Bucket 1 (100-199)
        model.add_value(350, vec![b'C'; 350], &mut rng);  // Bucket 3 (300-399)
        model.add_value(25000, vec![b'D'; 25000], &mut rng); // Catch-all bucket

        // Verify buckets contain the expected reads
        assert_eq!(model.batches[0].qualities.len(), 1);   // Bucket 0 has 1 read
        assert_eq!(model.batches[1].qualities.len(), 1);   // Bucket 1 has 1 read
        assert_eq!(model.batches[3].qualities.len(), 1);   // Bucket 3 has 1 read
        assert_eq!(model.batches.last().unwrap().qualities.len(), 1); // Catch-all has 1 read
    }

    #[test]
    fn test_reservoir_sampling_caps_bucket_size() {
        // Use small max capacity for faster test
        let mut model = QualityModel::new(Some(100), Some(20000), Some(10));
        let mut rng = StdRng::seed_from_u64(42);

        // Add 50 quality strings to the same bucket (bucket 0: 0-99)
        for _ in 0..50 {
            model.add_value(50, vec![b'X'; 50], &mut rng);
        }

        let batch = &model.batches[0];

        // Verify bucket is capped at max_capacity (10)
        assert_eq!(batch.qualities.len(), 10);

        // Verify total_seen tracks all 50 reads
        assert_eq!(batch.total_seen, 50);
    }

    #[test]
    fn test_truncation() {
        let mut model = QualityModel::new(None, None, None);
        let mut rng = StdRng::seed_from_u64(42);

        // Add a quality string of length 200
        model.add_value(200, vec![b'I'; 200], &mut rng);

        // Request length 150 - should be truncated from 200
        let sampled = model.sample(150, &mut rng);
        assert!(sampled.is_some());
        assert_eq!(sampled.unwrap().len(), 150);

        // Request length 200 - exact match
        let sampled = model.sample(200, &mut rng);
        assert!(sampled.is_some());
        assert_eq!(sampled.unwrap().len(), 200);

        // Request length 50 - should be truncated from 200
        let sampled = model.sample(50, &mut rng);
        assert!(sampled.is_some());
        assert_eq!(sampled.unwrap().len(), 50);
    }

    #[test]
    fn test_fallback_to_next_bucket() {
        let mut model = QualityModel::new(None, None, None);
        let mut rng = StdRng::seed_from_u64(42);

        // Only add quality string to bucket 5 (500-599)
        model.add_value(550, vec![b'K'; 550], &mut rng);

        // Request length 250 (bucket 2: 200-299) which is empty
        // Should fall back to bucket 5 and truncate to 250
        let sampled = model.sample(250, &mut rng);
        assert!(sampled.is_some());
        assert_eq!(sampled.unwrap().len(), 250);

        // Request length 400 (bucket 4: 400-499) which is also empty
        // Should fall back to bucket 5 and truncate to 400
        let sampled = model.sample(400, &mut rng);
        assert!(sampled.is_some());
        assert_eq!(sampled.unwrap().len(), 400);
    }

    #[test]
    fn test_empty_model_returns_none() {
        let model = QualityModel::new(None, None, None);
        let mut rng = StdRng::seed_from_u64(42);

        // Sampling from empty model should return None
        let sampled = model.sample(100, &mut rng);
        assert!(sampled.is_none());

        let sampled = model.sample(5000, &mut rng);
        assert!(sampled.is_none());
    }
}
