use rand::Rng;
use std::collections::HashMap;
use std::sync::LazyLock;

/// Mapping of Phred quality scores (0-93) to error probabilities.
///
/// Error probability is calculated as: 10^(-Q/10)
pub static QUALITY_MAPPING: LazyLock<HashMap<u8, f32>> = LazyLock::new(|| {
    (0..=93)
        .map(|q: u8| (q, 10.0_f32.powf(-(q as f32) / 10.0)))
        .collect()
});

/// Returns a random nucleotide different from the provided one.
///
/// # Arguments
/// * `nucleotide` - The nucleotide to exclude (as ASCII byte: b'A', b'C', b'G', or b'T')
/// * `rng` - Random number generator
///
/// # Returns
/// A random nucleotide byte from {A, C, G, T} excluding the input nucleotide
pub fn get_random_nucleotide<R: Rng>(nucleotide: u8, rng: &mut R) -> u8 {
    const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let idx = NUCLEOTIDES
        .iter()
        .position(|&n| n == nucleotide)
        .unwrap_or(0);
    let offset = rng.random_range(1..=3);

    NUCLEOTIDES[(idx + offset) % 4]
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_get_random_nucleotide() {
        let mut rng = StdRng::seed_from_u64(42);
        let result = get_random_nucleotide(b'A', &mut rng);
        assert_ne!(result, b'A');
        assert!(result == b'C' || result == b'G' || result == b'T');
    }

    #[test]
    fn test_quality_mapping() {
        // Q0 should be 1.0
        assert!((QUALITY_MAPPING[&0] - 1.0).abs() < 0.001);
        // Q10 should be approximately 0.1
        assert!((QUALITY_MAPPING[&10] - 0.1).abs() < 0.001);
        // Q20 should be approximately 0.01
        assert!((QUALITY_MAPPING[&20] - 0.01).abs() < 0.001);
    }
}
