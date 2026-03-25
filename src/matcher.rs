use memchr::memmem;

use crate::construct::Segment;

/// Result of matching a sequence against a library construct definition.
#[derive(Debug, Clone)]
pub struct MatchResult {
    /// Whether all segments matched and all rand regions are in range.
    pub matched: bool,
    /// Per-const segment match info (None = not found).
    pub const_matches: Vec<Option<ConstMatch>>,
    /// Per-rand segment match info (None = const before/after not found).
    pub rand_regions: Vec<Option<RandMatch>>,
    /// Leading gap: bases before first const segment.
    pub leading_gap: usize,
    /// Trailing gap: bases after last const segment.
    pub trailing_gap: usize,
    /// Failure reason, if any.
    pub failure: Option<FailureKind>,
}

#[derive(Debug, Clone)]
pub struct ConstMatch {
    pub position: usize,
    pub length: usize,
}

#[derive(Debug, Clone)]
pub struct RandMatch {
    pub start: usize,
    pub length: usize,
    pub in_range: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FailureKind {
    ConstNotFound(usize),
    RandOutOfRange(usize),
    LeadingGap(usize),
    TrailingGap(usize),
}

/// Precompiled matcher for a library construct definition.
pub struct ConstructMatcher {
    segments: Vec<Segment>,
    /// memmem::Finder for each const segment (indexed by const_index).
    finders: Vec<memmem::Finder<'static>>,
    /// Mapping: which indices in `segments` are Const.
    const_indices: Vec<usize>,
    /// Mapping: which indices in `segments` are Rand, with their position in const_indices.
    rand_specs: Vec<RandSpec>,
    strict_boundaries: bool,
}

struct RandSpec {
    segment_index: usize,
    min_len: usize,
    max_len: usize,
}

impl ConstructMatcher {
    pub fn new(segments: &[Segment], strict_boundaries: bool) -> Self {
        let mut finders = Vec::new();
        let mut const_indices = Vec::new();
        let mut rand_specs = Vec::new();

        for (i, seg) in segments.iter().enumerate() {
            match seg {
                Segment::Const { sequence } => {
                    // Leak the sequence to get a 'static Finder.
                    // This is intentional: ConstructMatcher lives for the program's lifetime.
                    let seq: &'static [u8] = Vec::leak(sequence.clone());
                    finders.push(memmem::Finder::new(seq));
                    const_indices.push(i);
                }
                Segment::Rand { min_len, max_len } => {
                    rand_specs.push(RandSpec {
                        segment_index: i,
                        min_len: *min_len,
                        max_len: *max_len,
                    });
                }
            }
        }

        Self {
            segments: segments.to_vec(),
            finders,
            const_indices,
            rand_specs,
            strict_boundaries,
        }
    }

    pub fn match_construct(&self, seq: &[u8]) -> MatchResult {
        let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        let seq_len = seq_upper.len();

        // 1. Find all const segments in order (left-to-right)
        let mut pos = 0;
        let mut anchors: Vec<Option<(usize, usize)>> = Vec::with_capacity(self.finders.len());

        for (finder_idx, finder) in self.finders.iter().enumerate() {
            match finder.find(&seq_upper[pos..]) {
                Some(offset) => {
                    let abs_start = pos + offset;
                    let const_seg_idx = self.const_indices[finder_idx];
                    let const_len = match &self.segments[const_seg_idx] {
                        Segment::Const { sequence } => sequence.len(),
                        _ => unreachable!(),
                    };
                    let abs_end = abs_start + const_len;
                    anchors.push(Some((abs_start, abs_end)));
                    pos = abs_end;
                }
                None => {
                    anchors.push(None);
                    // Build partial result
                    let const_matches = self.build_const_matches(&anchors);
                    return MatchResult {
                        matched: false,
                        const_matches,
                        rand_regions: Vec::new(),
                        leading_gap: 0,
                        trailing_gap: 0,
                        failure: Some(FailureKind::ConstNotFound(finder_idx)),
                    };
                }
            }
        }

        // All const segments found
        let const_matches = self.build_const_matches(&anchors);

        // 2. Extract gaps between anchors
        let anchor_positions: Vec<(usize, usize)> = anchors.iter().map(|a| a.unwrap()).collect();

        let leading_gap = anchor_positions[0].0;
        let trailing_gap = seq_len.saturating_sub(anchor_positions.last().unwrap().1);

        // 3. Map gaps to rand segments and validate
        let (rand_regions, failure) = self.validate_rand_regions(&anchor_positions, seq_len);

        // 4. Check boundaries
        let boundary_failure = if self.strict_boundaries {
            if leading_gap > 0 {
                Some(FailureKind::LeadingGap(leading_gap))
            } else if trailing_gap > 0 {
                Some(FailureKind::TrailingGap(trailing_gap))
            } else {
                None
            }
        } else {
            None
        };

        // Determine if construct starts/ends with rand
        let starts_with_rand = matches!(self.segments.first(), Some(Segment::Rand { .. }));
        let ends_with_rand = matches!(self.segments.last(), Some(Segment::Rand { .. }));

        // If construct starts with rand, leading gap is handled by rand validation
        // If construct ends with rand, trailing gap is handled by rand validation
        let effective_boundary_failure = if starts_with_rand && ends_with_rand {
            None
        } else if starts_with_rand {
            if trailing_gap > 0 && self.strict_boundaries {
                Some(FailureKind::TrailingGap(trailing_gap))
            } else {
                None
            }
        } else if ends_with_rand {
            if leading_gap > 0 && self.strict_boundaries {
                Some(FailureKind::LeadingGap(leading_gap))
            } else {
                None
            }
        } else {
            boundary_failure
        };

        let final_failure = failure.or(effective_boundary_failure);
        let matched = final_failure.is_none();

        MatchResult {
            matched,
            const_matches,
            rand_regions,
            leading_gap,
            trailing_gap,
            failure: final_failure,
        }
    }

    fn build_const_matches(&self, anchors: &[Option<(usize, usize)>]) -> Vec<Option<ConstMatch>> {
        anchors
            .iter()
            .map(|anchor| {
                anchor.map(|(start, end)| ConstMatch {
                    position: start,
                    length: end - start,
                })
            })
            .collect()
    }

    fn validate_rand_regions(
        &self,
        anchors: &[(usize, usize)],
        seq_len: usize,
    ) -> (Vec<Option<RandMatch>>, Option<FailureKind>) {
        let mut rand_regions = Vec::new();
        let mut failure = None;

        // Build all gaps: leading, inter-const, trailing
        let mut gaps = Vec::new();

        // Leading gap (before first const)
        gaps.push((0, anchors[0].0));

        // Inter-const gaps
        for i in 0..anchors.len() - 1 {
            gaps.push((anchors[i].1, anchors[i + 1].0));
        }

        // Trailing gap (after last const)
        gaps.push((anchors.last().unwrap().1, seq_len));

        // Map gaps to rand specs based on construct definition
        // Construct segments: [possibly rand] [const] [possibly rand] [const] ... [possibly rand]
        // gaps:               [leading]       [inter0]               [inter1] ... [trailing]
        //
        // If construct starts with Rand, leading gap maps to first rand spec
        // If construct ends with Rand, trailing gap maps to last rand spec

        let mut gap_idx = 0;
        let mut rand_spec_idx = 0;

        let starts_with_rand = matches!(self.segments.first(), Some(Segment::Rand { .. }));
        let ends_with_rand = matches!(self.segments.last(), Some(Segment::Rand { .. }));

        // Leading gap
        if starts_with_rand && rand_spec_idx < self.rand_specs.len() {
            let spec = &self.rand_specs[rand_spec_idx];
            let (start, end) = gaps[gap_idx];
            let len = end - start;
            let in_range = len >= spec.min_len && len <= spec.max_len;
            if !in_range && failure.is_none() {
                failure = Some(FailureKind::RandOutOfRange(rand_spec_idx));
            }
            rand_regions.push(Some(RandMatch {
                start,
                length: len,
                in_range,
            }));
            rand_spec_idx += 1;
        }
        gap_idx += 1;

        // Inter-const gaps
        while gap_idx < gaps.len() - 1 {
            if rand_spec_idx < self.rand_specs.len() {
                let spec = &self.rand_specs[rand_spec_idx];
                // Check if this rand spec belongs between these const segments
                let rand_seg_idx = spec.segment_index;
                let prev_const_seg_idx = self.const_indices[gap_idx - 1];
                let next_const_seg_idx = self.const_indices[gap_idx];
                if rand_seg_idx > prev_const_seg_idx && rand_seg_idx < next_const_seg_idx {
                    let (start, end) = gaps[gap_idx];
                    let len = end - start;
                    let in_range = len >= spec.min_len && len <= spec.max_len;
                    if !in_range && failure.is_none() {
                        failure = Some(FailureKind::RandOutOfRange(rand_spec_idx));
                    }
                    rand_regions.push(Some(RandMatch {
                        start,
                        length: len,
                        in_range,
                    }));
                    rand_spec_idx += 1;
                }
                // If no rand between these consts, the gap should be 0
                // (consecutive const segments with no rand in between)
            }
            gap_idx += 1;
        }

        // Trailing gap
        if ends_with_rand && rand_spec_idx < self.rand_specs.len() {
            let spec = &self.rand_specs[rand_spec_idx];
            let (start, end) = gaps.last().unwrap();
            let len = end - start;
            let in_range = len >= spec.min_len && len <= spec.max_len;
            if !in_range && failure.is_none() {
                failure = Some(FailureKind::RandOutOfRange(rand_spec_idx));
            }
            rand_regions.push(Some(RandMatch {
                start: *start,
                length: len,
                in_range,
            }));
        }

        (rand_regions, failure)
    }

    /// Number of const segments in this construct.
    pub fn const_count(&self) -> usize {
        self.finders.len()
    }

    /// Number of rand segments in this construct.
    pub fn rand_count(&self) -> usize {
        self.rand_specs.len()
    }

    /// Total number of segments.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Get segments.
    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }
}

/// Compute minimum Hamming distance of `pattern` against all windows of `seq`.
/// Used only for near-miss reporting on sequences that failed exact match.
pub fn near_miss_distance(seq: &[u8], pattern: &[u8]) -> usize {
    if seq.len() < pattern.len() {
        return usize::MAX;
    }
    seq.windows(pattern.len())
        .map(|window| {
            window
                .iter()
                .zip(pattern.iter())
                .filter(|(a, b)| a != b)
                .count()
        })
        .min()
        .unwrap_or(usize::MAX)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::construct::Segment;

    fn make_matcher(segments: &[Segment], strict: bool) -> ConstructMatcher {
        ConstructMatcher::new(segments, strict)
    }

    #[test]
    fn exact_match_simple() {
        let segments = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 5,
                max_len: 10,
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, true);
        let seq = b"ATCGNNNNNNGCTA";
        let result = matcher.match_construct(seq);
        assert!(result.matched);
        assert_eq!(result.const_matches.len(), 2);
        assert_eq!(result.const_matches[0].as_ref().unwrap().position, 0);
        assert_eq!(result.const_matches[1].as_ref().unwrap().position, 10);
        assert_eq!(result.rand_regions.len(), 1);
        assert_eq!(result.rand_regions[0].as_ref().unwrap().length, 6);
        assert!(result.rand_regions[0].as_ref().unwrap().in_range);
    }

    #[test]
    fn const_not_found() {
        let segments = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 5,
                max_len: 10,
            },
            Segment::Const {
                sequence: b"XXXX".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, true);
        let seq = b"ATCGNNNNNNGCTA";
        let result = matcher.match_construct(seq);
        assert!(!result.matched);
        assert_eq!(result.failure, Some(FailureKind::ConstNotFound(1)));
    }

    #[test]
    fn rand_out_of_range_too_short() {
        let segments = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 10,
                max_len: 20,
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, true);
        let seq = b"ATCGNNGCTA"; // only 2bp between consts
        let result = matcher.match_construct(seq);
        assert!(!result.matched);
        assert_eq!(result.failure, Some(FailureKind::RandOutOfRange(0)));
    }

    #[test]
    fn leading_gap_strict() {
        let segments = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 3,
                max_len: 5,
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, true);
        let seq = b"XXATCGNNNGCTA"; // 2bp leading gap
        let result = matcher.match_construct(seq);
        assert!(!result.matched);
        assert_eq!(result.leading_gap, 2);
        assert_eq!(result.failure, Some(FailureKind::LeadingGap(2)));
    }

    #[test]
    fn leading_gap_allowed() {
        let segments = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 3,
                max_len: 5,
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, false);
        let seq = b"XXATCGNNNGCTA";
        let result = matcher.match_construct(seq);
        assert!(result.matched);
        assert_eq!(result.leading_gap, 2);
    }

    #[test]
    fn leading_rand_absorbs_gap() {
        let segments = vec![
            Segment::Rand {
                min_len: 0,
                max_len: 5,
            },
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 3,
                max_len: 5,
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, true);
        let seq = b"XXATCGNNNGCTA"; // 2bp leading = within rand:0-5
        let result = matcher.match_construct(seq);
        assert!(result.matched);
    }

    #[test]
    fn consecutive_const_no_rand() {
        let segments = vec![
            Segment::Const {
                sequence: b"AAAA".to_vec(),
            },
            Segment::Const {
                sequence: b"TTTT".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, true);
        let seq = b"AAAATTTT";
        let result = matcher.match_construct(seq);
        assert!(result.matched);
        assert_eq!(result.const_matches.len(), 2);
    }

    #[test]
    fn case_insensitive() {
        let segments = vec![Segment::Const {
            sequence: b"ATCG".to_vec(),
        }];
        let matcher = make_matcher(&segments, false);
        let seq = b"atcg";
        let result = matcher.match_construct(seq);
        assert!(result.matched);
    }

    #[test]
    fn near_miss_distance_exact() {
        assert_eq!(near_miss_distance(b"ATCGATCG", b"ATCG"), 0);
    }

    #[test]
    fn near_miss_distance_one() {
        assert_eq!(near_miss_distance(b"AXCG", b"ATCG"), 1);
    }

    #[test]
    fn near_miss_distance_pattern_longer() {
        assert_eq!(near_miss_distance(b"AT", b"ATCG"), usize::MAX);
    }

    #[test]
    fn three_const_with_two_rand() {
        let segments = vec![
            Segment::Const {
                sequence: b"AAA".to_vec(),
            },
            Segment::Rand {
                min_len: 2,
                max_len: 4,
            },
            Segment::Const {
                sequence: b"TTT".to_vec(),
            },
            Segment::Rand {
                min_len: 1,
                max_len: 3,
            },
            Segment::Const {
                sequence: b"GGG".to_vec(),
            },
        ];
        let matcher = make_matcher(&segments, true);
        let seq = b"AAANNNTTTNNGGG";
        let result = matcher.match_construct(seq);
        assert!(result.matched, "failure: {:?}", result.failure);
        assert_eq!(result.rand_regions.len(), 2);
        assert_eq!(result.rand_regions[0].as_ref().unwrap().length, 3);
        assert_eq!(result.rand_regions[1].as_ref().unwrap().length, 2);
    }
}
