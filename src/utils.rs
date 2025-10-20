use rand::Rng;
use std::collections::HashMap;
use std::sync::LazyLock;

pub static QUALITY_MAPPING: LazyLock<HashMap<u8, f32>> = LazyLock::new(|| {
    (1..=50)
        .map(|q: u8| (q, 10.0_f32.powf(-(q as f32) / 10.0)))
        .collect()
});

pub fn get_random_nucleotide<R: Rng>(nucleotide: u8, mut rng: R) -> u8 {
    let nucleotides = [b'A', b'C', b'G', b'T'];
    let nucleotide_options: Vec<u8> = nucleotides
        .iter()
        .filter(|&&n| n != nucleotide)
        .copied()
        .collect();

    nucleotide_options[rng.random_range(0..nucleotide_options.len())]
}
