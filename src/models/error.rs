const SUBSTITUTION_DEFAULT_RATE: f64 = 0.7;
const INSERTION_DEFAULT_RATE: f64 = 0.1;
const DELETION_DEFAULT_RATE: f64 = 0.2;

pub struct ErrorModel {
    pub substitution_rate: f64,
    pub insertion_rate: f64,
    pub deletion_rate: f64,
}

impl ErrorModel {
    pub fn new(
        substitution_rate: Option<f64>,
        insertion_rate: Option<f64>,
        deletion_rate: Option<f64>,
    ) -> Self {
        ErrorModel {
            substitution_rate: substitution_rate.unwrap_or(SUBSTITUTION_DEFAULT_RATE),
            insertion_rate: insertion_rate.unwrap_or(INSERTION_DEFAULT_RATE),
            deletion_rate: deletion_rate.unwrap_or(DELETION_DEFAULT_RATE),
        }
    }
}
