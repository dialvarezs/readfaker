use crate::io::fastq::FastqReader;
use crate::models::{LengthModel, QualityModel};
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::path::Path;
use std::sync::LazyLock;

/// Mapping of Phred quality scores (0-93) to error probabilities.
///
/// Error probability is calculated as: 10^(-Q/10)
pub static QUALITY_MAPPING: LazyLock<[f32; 94]> = LazyLock::new(|| {
    let mut mapping = [0.0_f32; 94];
    for (q, slot) in mapping.iter_mut().enumerate() {
        *slot = 10.0_f32.powf(-(q as f32) / 10.0);
    }
    mapping
});

/// Loads length and quality models from an existing FASTQ file.
///
/// Reads all records from the input file and builds empirical models
/// for read lengths and quality scores.
///
/// # Arguments
/// * `fastq_path` - Path to the FASTQ file to analyze
///
/// # Returns
/// Tuple of (LengthModel, QualityModel) built from the input file
pub fn load_models(
    fastq_path: &Path,
    seed: Option<u64>,
) -> anyhow::Result<(LengthModel, QualityModel)> {
    let mut length_model = LengthModel::new();
    let mut quality_model = QualityModel::new(None, None, None);

    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_rng(&mut rand::rng()),
    };

    let reader = FastqReader::from_path(fastq_path)?;

    for record in reader {
        let record = record?;
        length_model.add_value(record.len());
        quality_model.add_value(record.len(), record.quality, &mut rng);
    }

    Ok((length_model, quality_model))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_mapping() {
        // Q0 should be 1.0
        assert!((QUALITY_MAPPING[0] - 1.0).abs() < 0.001);
        // Q10 should be approximately 0.1
        assert!((QUALITY_MAPPING[10] - 0.1).abs() < 0.001);
        // Q20 should be approximately 0.01
        assert!((QUALITY_MAPPING[20] - 0.01).abs() < 0.001);
    }
}
