use crate::io::bam::BamReader;
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

    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_rng(&mut rand::rng()),
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
            let length = record.sequence().len();
            // Convert raw Phred scores (0-93) to Phred+33 ASCII encoding
            let quality: Vec<u8> = record
                .quality_scores()
                .as_ref()
                .iter()
                .map(|&q| q.saturating_add(33))
                .collect();
            length_model.add_value(length);
            quality_model.add_value(length, quality, &mut rng);
        }
    }
    else if base_name.ends_with(".fastq") || base_name.ends_with(".fq") {
        let reader = FastqReader::from_path(input_path)?;
        for record in reader {
            let record = record?;
            let length = record.sequence().len();
            let quality = record.quality_scores().to_vec();
            length_model.add_value(length);
            quality_model.add_value(length, quality, &mut rng);
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
    fn test_quality_mapping() {
        // Q0 should be 1.0
        assert!((QUALITY_MAPPING[0] - 1.0).abs() < 0.001);
        // Q10 should be approximately 0.1
        assert!((QUALITY_MAPPING[10] - 0.1).abs() < 0.001);
        // Q20 should be approximately 0.01
        assert!((QUALITY_MAPPING[20] - 0.01).abs() < 0.001);
    }
}
