use anyhow::Result;
use std::path::PathBuf;
use crate::distributions::{LengthDistribution, QualityDistribution};
use crate::io::fastq::FastqReader;

pub fn generate_reads() {

}

pub fn load_distributions(fastq_path: PathBuf) -> Result<(LengthDistribution, QualityDistribution)> {
    let mut length_distribution = LengthDistribution::new();
    let mut quality_distribution = QualityDistribution::new();

    let reader = FastqReader::from_path(&fastq_path)?;

    for record in reader {
        let record = record?;
        length_distribution.add_value(record.len());
        quality_distribution.add_value(record.len(), record.quality);
    }

    Ok((length_distribution, quality_distribution))
}