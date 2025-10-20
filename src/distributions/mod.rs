//! Statistical distributions for read length and quality scores.
//!
//! These distributions are built from empirical data and used to generate
//! realistic synthetic reads.

pub mod length;
pub mod quality;

pub use length::LengthDistribution;
pub use quality::QualityDistribution;
