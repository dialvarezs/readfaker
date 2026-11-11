use crate::io::bam::BamReader;
use crate::io::fastq::FastqReader;
use crate::models::{LengthModel, QualityModel};
use rand::SeedableRng;
use rand::rngs::SmallRng;
use std::path::Path;
use std::sync::LazyLock;

/// Offset of Phred+33 ASCII quality encoding.
pub const PHRED_OFFSET: u8 = 33;

/// Per-base error thresholds indexed by Phred+33 ASCII byte (Q clamped to 0-93).
///
/// Error probability 10^(-Q/10) scaled to the `u32` range, so an error occurs when a
/// uniform random `u32` is <= the threshold. Avoids float math in the per-base hot loop.
pub static ERROR_THRESHOLDS: LazyLock<[u32; 256]> = LazyLock::new(|| {
    std::array::from_fn(|ascii| {
        let q = (ascii as u8).saturating_sub(PHRED_OFFSET).min(93);
        (10.0_f64.powf(-f64::from(q) / 10.0) * f64::from(u32::MAX)) as u32
    })
});

/// Creates a seeded RNG, or one seeded from system entropy if `seed` is None.
pub fn make_rng(seed: Option<u64>) -> SmallRng {
    match seed {
        Some(s) => SmallRng::seed_from_u64(s),
        None => SmallRng::from_rng(&mut rand::rng()),
    }
}

/// Loads length and quality models from an existing FASTQ or BAM file.
///
/// Automatically detects the file format based on the extension (.fastq, .fq, .bam).
/// Reads all records from the input file and builds empirical models
/// for read lengths and quality scores.
///
/// # Arguments
/// * `input_path` - Path to the FASTQ or BAM file to analyze
/// * `seed` - Optional random seed for reproducibility (uses system entropy if None)
///
/// # Returns
/// Tuple of (LengthModel, QualityModel) built from the input file
pub fn load_models(
    input_path: &Path,
    seed: Option<u64>,
) -> anyhow::Result<(LengthModel, QualityModel)> {
    let mut length_model = LengthModel::new();
    let mut quality_model = QualityModel::new(None, None, None);

    let mut rng = make_rng(seed);
    let mut add = |length: usize, quality: Vec<u8>| {
        length_model.add_value(length);
        quality_model.add_value(length, quality, &mut rng);
    };

    // Detect file format based on file name (handles compound extensions like .fastq.gz)
    let file_name = input_path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");

    let file_name_lower = file_name.to_lowercase();

    // Strip compression extensions to get base extension
    let base_name = file_name_lower
        .strip_suffix(".gz")
        .or_else(|| file_name_lower.strip_suffix(".bgz"))
        .unwrap_or(&file_name_lower);

    if base_name.ends_with(".bam") {
        let reader = BamReader::from_path(input_path)?;
        for record in reader {
            let record = record?;
            // Convert raw Phred scores (0-93) to Phred+33 ASCII encoding
            let quality = record
                .quality_scores()
                .as_ref()
                .iter()
                .map(|&q| q.saturating_add(PHRED_OFFSET))
                .collect();
            add(record.sequence().len(), quality);
        }
    } else if base_name.ends_with(".fastq") || base_name.ends_with(".fq") {
        let reader = FastqReader::from_path(input_path)?;
        for record in reader {
            let record = record?;
            add(record.sequence().len(), record.quality_scores().to_vec());
        }
    } else {
        anyhow::bail!(
            "Unsupported input file format. Expected .fastq, .fq, or .bam (optionally compressed with .gz, .bgz)"
        );
    }

    Ok((length_model, quality_model))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_thresholds() {
        let p = |ascii: u8| f64::from(ERROR_THRESHOLDS[usize::from(ascii)]) / f64::from(u32::MAX);
        // Q0 always errors, Q10 ~ 0.1, Q20 ~ 0.01 (ASCII = Q + 33)
        assert_eq!(ERROR_THRESHOLDS[usize::from(b'!')], u32::MAX);
        assert!((p(b'+') - 0.1).abs() < 0.001);
        assert!((p(b'5') - 0.01).abs() < 0.001);
        // Bytes below the offset clamp to Q0
        assert_eq!(ERROR_THRESHOLDS[0], u32::MAX);
    }
}
