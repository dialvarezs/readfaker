use rand::Rng;
use std::collections::HashMap;
use std::sync::LazyLock;

/// Mapping of Phred quality scores (1-50) to error probabilities.
///
/// Error probability is calculated as: 10^(-Q/10)
pub static QUALITY_MAPPING: LazyLock<HashMap<u8, f32>> = LazyLock::new(|| {
    (1..=50)
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
pub fn get_random_nucleotide<R: Rng>(nucleotide: u8, mut rng: R) -> u8 {
    let nucleotides = [b'A', b'C', b'G', b'T'];
    let nucleotide_options: Vec<u8> = nucleotides
        .iter()
        .filter(|&&n| n != nucleotide)
        .copied()
        .collect();

    nucleotide_options[rng.random_range(0..nucleotide_options.len())]
}
