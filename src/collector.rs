use ahash::AHashMap;
use seqtable::DualSeqCounts;

use crate::construct::Segment;
use crate::matcher::{MatchResult, near_miss_distance};

/// Sequence frequency counting, backed by seqtable's DualSeqCounts.
pub struct CountCollector {
    counts: DualSeqCounts,
    total_reads: u64,
}

pub struct CountEntry {
    pub sequence: Vec<u8>,
    pub count: u64,
    pub rpm: f64,
    pub rank: u64,
}

impl Default for CountCollector {
    fn default() -> Self {
        Self::new()
    }
}

impl CountCollector {
    pub fn new() -> Self {
        Self {
            counts: DualSeqCounts::new(),
            total_reads: 0,
        }
    }

    /// Create from pre-counted DualSeqCounts (for count-only mode).
    pub fn from_counts(counts: DualSeqCounts, total_reads: u64) -> Self {
        Self {
            counts,
            total_reads,
        }
    }

    pub fn process_record(&mut self, seq: &[u8]) {
        self.counts.insert(seq);
        self.total_reads += 1;
    }

    pub fn total_reads(&self) -> u64 {
        self.total_reads
    }

    pub fn unique_count(&self) -> usize {
        self.counts.len()
    }

    pub fn finalize(self) -> Vec<CountEntry> {
        let total = self.total_reads as f64;
        let mut entries: Vec<CountEntry> = self
            .counts
            .into_iter()
            .map(|(sequence, count)| CountEntry {
                sequence,
                count,
                rpm: if total > 0.0 {
                    (count as f64 / total) * 1_000_000.0
                } else {
                    0.0
                },
                rank: 0,
            })
            .collect();

        entries.sort_unstable_by(|a, b| b.count.cmp(&a.count));

        // Competition ranking (1-2-2-4)
        let mut rank = 1u64;
        let mut i = 0;
        while i < entries.len() {
            let current_count = entries[i].count;
            let mut j = i;
            while j < entries.len() && entries[j].count == current_count {
                entries[j].rank = rank;
                j += 1;
            }
            rank = j as u64 + 1;
            i = j;
        }

        entries
    }
}

/// Per-const segment match statistics + near-miss analysis.
pub struct ConstructCollector {
    total: u64,
    valid: u64,
    const_stats: Vec<ConstSegmentStats>,
    rand_stats: Vec<RandSegmentStats>,
    const_sequences: Vec<Vec<u8>>,
    near_miss_depth: usize,
}

pub struct ConstSegmentStats {
    pub found: u64,
    pub not_found: u64,
    pub near_miss_dist: [u64; 4], // 1mm, 2mm, 3mm, 4+mm
}

pub struct RandSegmentStats {
    pub length_dist: AHashMap<usize, u64>,
    pub in_range: u64,
    pub out_of_range: u64,
}

impl ConstructCollector {
    pub fn new(segments: &[Segment], near_miss_depth: usize) -> Self {
        let const_count = segments
            .iter()
            .filter(|s| matches!(s, Segment::Const { .. }))
            .count();
        let rand_count = segments
            .iter()
            .filter(|s| matches!(s, Segment::Rand { .. }))
            .count();

        let const_sequences: Vec<Vec<u8>> = segments
            .iter()
            .filter_map(|s| match s {
                Segment::Const { sequence } => Some(sequence.clone()),
                _ => None,
            })
            .collect();

        Self {
            total: 0,
            valid: 0,
            const_stats: (0..const_count)
                .map(|_| ConstSegmentStats {
                    found: 0,
                    not_found: 0,
                    near_miss_dist: [0; 4],
                })
                .collect(),
            rand_stats: (0..rand_count)
                .map(|_| RandSegmentStats {
                    length_dist: AHashMap::new(),
                    in_range: 0,
                    out_of_range: 0,
                })
                .collect(),
            const_sequences,
            near_miss_depth,
        }
    }

    pub fn process_record(&mut self, seq: &[u8], result: &MatchResult) {
        self.total += 1;
        if result.matched {
            self.valid += 1;
        }

        // Per-const stats
        for (i, cm) in result.const_matches.iter().enumerate() {
            if i >= self.const_stats.len() {
                break;
            }
            match cm {
                Some(_) => {
                    self.const_stats[i].found += 1;
                }
                None => {
                    self.const_stats[i].not_found += 1;
                    if self.near_miss_depth > 0 && i < self.const_sequences.len() {
                        let seq_upper: Vec<u8> =
                            seq.iter().map(|b| b.to_ascii_uppercase()).collect();
                        let dist = near_miss_distance(&seq_upper, &self.const_sequences[i]);
                        if dist <= self.near_miss_depth {
                            let bucket = (dist - 1).min(3);
                            self.const_stats[i].near_miss_dist[bucket] += 1;
                        }
                    }
                }
            }
        }

        // Per-rand stats
        for (i, rm) in result.rand_regions.iter().enumerate() {
            if i >= self.rand_stats.len() {
                break;
            }
            if let Some(rm) = rm {
                *self.rand_stats[i].length_dist.entry(rm.length).or_insert(0) += 1;
                if rm.in_range {
                    self.rand_stats[i].in_range += 1;
                } else {
                    self.rand_stats[i].out_of_range += 1;
                }
            }
        }
    }

    pub fn total(&self) -> u64 {
        self.total
    }

    pub fn valid(&self) -> u64 {
        self.valid
    }

    pub fn const_stats(&self) -> &[ConstSegmentStats] {
        &self.const_stats
    }

    pub fn rand_stats(&self) -> &[RandSegmentStats] {
        &self.rand_stats
    }
}

/// Per-position base composition for each rand segment.
pub struct CompositionCollector {
    /// For each rand segment: per-position [A, T, G, C] counts.
    /// Outer Vec = rand segments, inner Vec = positions.
    per_rand: Vec<Vec<[u64; 4]>>,
    rand_count: usize,
}

impl CompositionCollector {
    pub fn new(segments: &[Segment]) -> Self {
        let rand_count = segments
            .iter()
            .filter(|s| matches!(s, Segment::Rand { .. }))
            .count();
        Self {
            per_rand: vec![Vec::new(); rand_count],
            rand_count,
        }
    }

    pub fn process_record(&mut self, seq: &[u8], result: &MatchResult) {
        if !result.matched {
            return;
        }
        let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();

        for (i, rm) in result.rand_regions.iter().enumerate() {
            if i >= self.rand_count {
                break;
            }
            if let Some(rm) = rm {
                let region = &seq_upper[rm.start..rm.start + rm.length];
                // Extend position array if needed
                while self.per_rand[i].len() < region.len() {
                    self.per_rand[i].push([0; 4]);
                }
                for (pos, &base) in region.iter().enumerate() {
                    let idx = match base {
                        b'A' => 0,
                        b'T' | b'U' => 1,
                        b'G' => 2,
                        b'C' => 3,
                        _ => continue,
                    };
                    self.per_rand[i][pos][idx] += 1;
                }
            }
        }
    }

    pub fn per_rand(&self) -> &[Vec<[u64; 4]>] {
        &self.per_rand
    }
}

/// Bitmask-based combination analysis for all segment conditions.
pub struct CombinationCollector {
    combinations: AHashMap<u32, u64>,
    num_const: usize,
    num_rand: usize,
}

impl CombinationCollector {
    pub fn new(segments: &[Segment]) -> Self {
        let num_const = segments
            .iter()
            .filter(|s| matches!(s, Segment::Const { .. }))
            .count();
        let num_rand = segments
            .iter()
            .filter(|s| matches!(s, Segment::Rand { .. }))
            .count();
        Self {
            combinations: AHashMap::new(),
            num_const,
            num_rand,
        }
    }

    pub fn process_record(&mut self, result: &MatchResult) {
        let mut mask: u32 = 0;

        for (i, cm) in result.const_matches.iter().enumerate() {
            if cm.is_some() {
                mask |= 1 << i;
            }
        }

        for (i, rm) in result.rand_regions.iter().enumerate() {
            if rm.as_ref().is_some_and(|r| r.in_range) {
                mask |= 1 << (self.num_const + i);
            }
        }

        *self.combinations.entry(mask).or_insert(0) += 1;
    }

    /// All-pass mask (every segment condition met).
    pub fn all_pass_mask(&self) -> u32 {
        (1u32 << (self.num_const + self.num_rand)) - 1
    }

    pub fn combinations(&self) -> &AHashMap<u32, u64> {
        &self.combinations
    }

    pub fn num_const(&self) -> usize {
        self.num_const
    }

    pub fn num_rand(&self) -> usize {
        self.num_rand
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::construct::Segment;
    use crate::matcher::ConstructMatcher;

    fn make_result(segments: &[Segment], seq: &[u8]) -> MatchResult {
        let matcher = ConstructMatcher::new(segments, true);
        matcher.match_construct(seq)
    }

    #[test]
    fn count_collector_basic() {
        let mut col = CountCollector::new();
        col.process_record(b"AATT");
        col.process_record(b"AATT");
        col.process_record(b"AACC");
        assert_eq!(col.total_reads(), 3);
        assert_eq!(col.unique_count(), 2);

        let entries = col.finalize();
        assert_eq!(entries[0].count, 2);
        assert_eq!(entries[0].rank, 1);
        assert_eq!(entries[1].count, 1);
        assert_eq!(entries[1].rank, 2);
    }

    #[test]
    fn construct_collector_found_not_found() {
        let segments = vec![
            Segment::Const {
                sequence: b"ATCG".to_vec(),
            },
            Segment::Rand {
                min_len: 2,
                max_len: 5,
            },
            Segment::Const {
                sequence: b"GCTA".to_vec(),
            },
        ];
        let mut col = ConstructCollector::new(&segments, 3);

        let r1 = make_result(&segments, b"ATCGNNGCTA");
        col.process_record(b"ATCGNNGCTA", &r1);
        let r2 = make_result(&segments, b"ATCGNNXXXXX");
        col.process_record(b"ATCGNNXXXXX", &r2);

        assert_eq!(col.total(), 2);
        assert_eq!(col.valid(), 1);
        assert_eq!(col.const_stats()[0].found, 2);
        assert_eq!(col.const_stats()[1].found, 1);
        assert_eq!(col.const_stats()[1].not_found, 1);
    }

    #[test]
    fn composition_collector_basic() {
        let segments = vec![
            Segment::Const {
                sequence: b"AA".to_vec(),
            },
            Segment::Rand {
                min_len: 3,
                max_len: 3,
            },
            Segment::Const {
                sequence: b"TT".to_vec(),
            },
        ];
        let mut col = CompositionCollector::new(&segments);

        for seq in [b"AAATCTT" as &[u8], b"AAATCTT"] {
            let r = make_result(&segments, seq);
            col.process_record(seq, &r);
        }

        let per_rand = col.per_rand();
        assert_eq!(per_rand.len(), 1);
        assert_eq!(per_rand[0].len(), 3);
        assert_eq!(per_rand[0][0][0], 2); // A
        assert_eq!(per_rand[0][1][1], 2); // T
        assert_eq!(per_rand[0][2][3], 2); // C
    }

    #[test]
    fn combination_collector_bitmask() {
        let segments = vec![
            Segment::Const {
                sequence: b"AA".to_vec(),
            },
            Segment::Rand {
                min_len: 2,
                max_len: 4,
            },
            Segment::Const {
                sequence: b"TT".to_vec(),
            },
        ];
        let mut col = CombinationCollector::new(&segments);

        for seq in [b"AANNTT" as &[u8], b"AANNXXXX"] {
            let r = make_result(&segments, seq);
            col.process_record(&r);
        }

        assert_eq!(*col.combinations().get(&0b111).unwrap_or(&0), 1);
        assert_eq!(*col.combinations().get(&0b001).unwrap_or(&0), 1);
        assert_eq!(col.all_pass_mask(), 0b111);
    }
}
