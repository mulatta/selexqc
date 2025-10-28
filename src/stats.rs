use serde::{Deserialize, Serialize, Serializer};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FailureReason {
    MissingConstant,
    IncorrectTotalLength,
    IncorrectStructure,
    LowQuality,
}

impl std::fmt::Display for FailureReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FailureReason::MissingConstant => write!(f, "Missing constant region"),
            FailureReason::IncorrectTotalLength => write!(f, "Incorrect total length"),
            FailureReason::IncorrectStructure => {
                write!(f, "Incorrect structure (upstream/downstream)")
            }
            FailureReason::LowQuality => write!(f, "Low quality score"),
        }
    }
}

fn serialize_pair_distribution<S>(
    map: &HashMap<(usize, usize), usize>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    use serde::ser::SerializeMap;
    let mut ser_map = serializer.serialize_map(Some(map.len()))?;
    for ((up, down), count) in map {
        let key = format!("({},{})", up, down);
        ser_map.serialize_entry(&key, count)?;
    }
    ser_map.end()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationConfig {
    pub constant_region: String,
    #[serde(skip)]
    pub validation_mode: crate::validate::ValidationMode,
    pub min_length: Option<usize>,
    pub max_length: Option<usize>,
    pub upstream_length: Option<usize>,
    pub upstream_tolerance: Option<usize>,
    pub downstream_length: Option<usize>,
    pub downstream_tolerance: Option<usize>,
    pub min_quality: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Stats {
    pub total_sequences: usize,
    pub valid_sequences: usize,
    pub invalid_sequences: usize,
    pub has_constant_region: usize,
    pub missing_constant_region: usize,
    pub correct_total_length: usize,
    pub incorrect_total_length: usize,
    pub correct_upstream: usize,
    pub correct_downstream: usize,
    pub correct_structure_pair: usize,
    pub correct_quality: usize,
    pub length_distribution: HashMap<usize, usize>,
    pub upstream_distribution: HashMap<usize, usize>,
    pub downstream_distribution: HashMap<usize, usize>,
    #[serde(serialize_with = "serialize_pair_distribution")]
    pub structure_pair_distribution: HashMap<(usize, usize), usize>,
    pub failure_statistics: HashMap<FailureReason, usize>,
    pub config: ValidationConfig,
}

impl Stats {
    pub fn valid_pct(&self) -> f64 {
        percentage(self.valid_sequences, self.total_sequences)
    }

    pub fn constant_pct(&self) -> f64 {
        percentage(self.has_constant_region, self.total_sequences)
    }

    pub fn correct_length_pct(&self) -> f64 {
        percentage(self.correct_total_length, self.total_sequences)
    }

    pub fn correct_upstream_pct(&self) -> f64 {
        percentage(self.correct_upstream, self.has_constant_region)
    }

    pub fn correct_downstream_pct(&self) -> f64 {
        percentage(self.correct_downstream, self.has_constant_region)
    }

    pub fn correct_structure_pair_pct(&self) -> f64 {
        percentage(self.correct_structure_pair, self.has_constant_region)
    }

    pub fn correct_quality_pct(&self) -> f64 {
        percentage(self.correct_quality, self.total_sequences)
    }

    pub fn filtered_pct(&self) -> f64 {
        100.0 - self.valid_pct()
    }
}

fn percentage(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        (numerator as f64 / denominator as f64) * 100.0
    }
}
