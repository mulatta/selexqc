use anyhow::{Result, bail};

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Segment {
    Const { sequence: Vec<u8> },
    Rand { min_len: usize, max_len: usize },
}

/// Parse a segment from CLI string: "const:ATCG" or "rand:35-45"
pub fn parse_segment(s: &str) -> Result<Segment> {
    let (kind, value) = s
        .split_once(':')
        .ok_or_else(|| anyhow::anyhow!("Expected 'const:SEQ' or 'rand:MIN-MAX', got '{}'", s))?;

    match kind {
        "const" => {
            let seq = value.to_uppercase();
            validate_nucleotide_sequence(&seq)?;
            Ok(Segment::Const {
                sequence: seq.into_bytes(),
            })
        }
        "rand" => {
            let (min_s, max_s) = value
                .split_once('-')
                .ok_or_else(|| anyhow::anyhow!("Expected 'rand:MIN-MAX', got 'rand:{}'", value))?;
            let min_len: usize = min_s.parse()?;
            let max_len: usize = max_s.parse()?;
            if min_len > max_len {
                bail!("rand min ({}) > max ({})", min_len, max_len);
            }
            Ok(Segment::Rand { min_len, max_len })
        }
        _ => bail!("Unknown segment type '{}'. Use 'const' or 'rand'", kind),
    }
}

fn validate_nucleotide_sequence(seq: &str) -> Result<()> {
    if seq.is_empty() {
        bail!("const sequence cannot be empty");
    }
    for (i, c) in seq.chars().enumerate() {
        if !matches!(
            c,
            'A' | 'T'
                | 'C'
                | 'G'
                | 'U'
                | 'R'
                | 'Y'
                | 'M'
                | 'K'
                | 'S'
                | 'W'
                | 'B'
                | 'D'
                | 'H'
                | 'V'
                | 'N'
        ) {
            bail!("Invalid nucleotide '{}' at position {} in '{}'", c, i, seq);
        }
    }
    Ok(())
}

/// Validate segment ordering rules:
/// - At least one const segment
/// - No consecutive rand segments
pub fn validate_segments(segments: &[Segment]) -> Result<()> {
    if segments.is_empty() {
        bail!("At least one segment is required");
    }

    let has_const = segments.iter().any(|s| matches!(s, Segment::Const { .. }));
    if !has_const {
        bail!("At least one const segment is required");
    }

    for window in segments.windows(2) {
        if matches!(
            (&window[0], &window[1]),
            (Segment::Rand { .. }, Segment::Rand { .. })
        ) {
            bail!("Consecutive rand segments are not allowed (ambiguous boundaries)");
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_const_segment() {
        let seg = parse_segment("const:ATCGATCG").unwrap();
        assert_eq!(
            seg,
            Segment::Const {
                sequence: b"ATCGATCG".to_vec()
            }
        );
    }

    #[test]
    fn parse_const_lowercase() {
        let seg = parse_segment("const:atcg").unwrap();
        assert_eq!(
            seg,
            Segment::Const {
                sequence: b"ATCG".to_vec()
            }
        );
    }

    #[test]
    fn parse_rand_segment() {
        let seg = parse_segment("rand:35-45").unwrap();
        assert_eq!(
            seg,
            Segment::Rand {
                min_len: 35,
                max_len: 45
            }
        );
    }

    #[test]
    fn parse_rand_invalid_range() {
        assert!(parse_segment("rand:45-35").is_err());
    }

    #[test]
    fn parse_invalid_type() {
        assert!(parse_segment("foo:bar").is_err());
    }

    #[test]
    fn parse_missing_colon() {
        assert!(parse_segment("constATCG").is_err());
    }

    #[test]
    fn parse_empty_const() {
        assert!(parse_segment("const:").is_err());
    }

    #[test]
    fn parse_invalid_nucleotide() {
        assert!(parse_segment("const:ATCGX").is_err());
    }

    #[test]
    fn parse_iupac_ambiguity_codes() {
        let seg = parse_segment("const:RYMKSWBDHVN").unwrap();
        assert!(matches!(seg, Segment::Const { .. }));
    }

    #[test]
    fn validate_empty() {
        assert!(validate_segments(&[]).is_err());
    }

    #[test]
    fn validate_no_const() {
        let segs = vec![Segment::Rand {
            min_len: 10,
            max_len: 20,
        }];
        assert!(validate_segments(&segs).is_err());
    }

    #[test]
    fn validate_consecutive_rand() {
        let segs = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 10,
                max_len: 20,
            },
            Segment::Rand {
                min_len: 5,
                max_len: 10,
            },
        ];
        assert!(validate_segments(&segs).is_err());
    }

    #[test]
    fn validate_valid_construct() {
        let segs = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 35,
                max_len: 45,
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        assert!(validate_segments(&segs).is_ok());
    }

    #[test]
    fn validate_consecutive_const_ok() {
        let segs = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        assert!(validate_segments(&segs).is_ok());
    }
}
