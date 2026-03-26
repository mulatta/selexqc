use std::path::PathBuf;

use anyhow::{Context, Result};
use rayon::prelude::*;

use crate::io;

pub struct ClusterConfig {
    pub input: PathBuf,
    pub output_prefix: String,
    pub max_dist: u32,
    pub count_filter: Option<u64>,
    pub max_clusters: Option<usize>,
    pub quiet: bool,
}

pub struct ClusterResult {
    pub cluster_id: u32,
    pub members: Vec<ClusterMember>,
    pub total_reads: u64,
}

pub struct ClusterMember {
    pub sequence: String,
    pub count: u64,
    pub rpm: f64,
    pub original_rank: u64,
    pub edit_distance: u32,
    pub cluster_rank: u32,
}

pub fn run_cluster(config: &ClusterConfig) -> Result<()> {
    let path = config.input.to_str().unwrap();
    if !config.quiet {
        eprintln!("  Reading: {}", path);
    }

    let entries = io::read_count_parquet(path)?;

    // Apply count filter
    let entries: Vec<_> = match config.count_filter {
        Some(min_count) => entries.into_iter().filter(|e| e.count >= min_count).collect(),
        None => entries,
    };

    if entries.is_empty() {
        anyhow::bail!("No sequences to cluster after filtering");
    }

    if !config.quiet {
        eprintln!(
            "  Clustering {} sequences (max edit distance: {})",
            entries.len(),
            config.max_dist
        );
    }

    let clusters = greedy_cluster(&entries, config.max_dist, config.max_clusters, config.quiet);

    if !config.quiet {
        let total_clustered: usize = clusters.iter().map(|c| c.members.len()).sum();
        eprintln!(
            "  {} clusters from {} sequences",
            clusters.len(),
            total_clustered
        );
    }

    // Write outputs
    let fasta_path = format!("{}.clustered.fasta", config.output_prefix);
    write_clustered_fasta(&fasta_path, &clusters)?;

    let summary_path = format!("{}.clusters.tsv", config.output_prefix);
    write_cluster_summary(&summary_path, &clusters)?;

    Ok(())
}

fn greedy_cluster(
    entries: &[io::ParquetEntry],
    max_dist: u32,
    max_clusters: Option<usize>,
    quiet: bool,
) -> Vec<ClusterResult> {
    // Pool of indices into entries (sorted by count desc, already from Parquet)
    let mut pool: Vec<usize> = (0..entries.len()).collect();
    let mut clusters = Vec::new();
    let mut cluster_id = 0u32;

    let progress = if !quiet && entries.len() > 1000 {
        let pb = indicatif::ProgressBar::new(entries.len() as u64);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("  Clustering [{bar:40}] {pos}/{len} ({percent}%)")
                .unwrap()
                .progress_chars("##-"),
        );
        Some(pb)
    } else {
        None
    };

    while !pool.is_empty() {
        if let Some(max) = max_clusters
            && clusters.len() >= max {
                break;
            }

        // Seed = first in pool (highest count among unassigned)
        let seed_idx = pool[0];
        let seed = &entries[seed_idx];
        let seed_bytes = seed.sequence.as_bytes();
        cluster_id += 1;

        // Parallel: compute distance for all remaining pool entries
        let pool_rest = &pool[1..];
        let matched: Vec<(usize, u32)> = pool_rest
            .par_iter()
            .filter_map(|&idx| {
                let candidate = &entries[idx];
                levenshtein_bounded(seed_bytes, candidate.sequence.as_bytes(), max_dist)
                    .map(|d| (idx, d))
            })
            .collect();

        // Build cluster
        let mut members = Vec::with_capacity(matched.len() + 1);
        members.push(ClusterMember {
            sequence: seed.sequence.clone(),
            count: seed.count,
            rpm: seed.rpm,
            original_rank: seed.rank,
            edit_distance: 0,
            cluster_rank: 1,
        });

        // Collect matched indices for removal from pool
        let mut matched_set: Vec<bool> = vec![false; entries.len()];
        matched_set[seed_idx] = true;

        for &(idx, dist) in &matched {
            matched_set[idx] = true;
            members.push(ClusterMember {
                sequence: entries[idx].sequence.clone(),
                count: entries[idx].count,
                rpm: entries[idx].rpm,
                original_rank: entries[idx].rank,
                edit_distance: dist,
                cluster_rank: 0, // assigned below
            });
        }

        // Assign cluster_rank by count descending
        members.sort_by(|a, b| b.count.cmp(&a.count));
        for (i, m) in members.iter_mut().enumerate() {
            m.cluster_rank = i as u32 + 1;
        }

        let total_reads = members.iter().map(|m| m.count).sum();

        if let Some(ref pb) = progress {
            pb.inc(members.len() as u64);
        }

        clusters.push(ClusterResult {
            cluster_id,
            members,
            total_reads,
        });

        // Remove matched from pool
        pool.retain(|idx| !matched_set[*idx]);
    }

    if let Some(pb) = progress {
        pb.finish_and_clear();
    }

    clusters
}

/// Bounded banded Levenshtein distance.
/// Returns None if distance > max_dist (early termination).
pub fn levenshtein_bounded(a: &[u8], b: &[u8], max_dist: u32) -> Option<u32> {
    let (m, n) = (a.len(), b.len());

    // Length difference alone exceeds threshold
    let len_diff = (m as i64 - n as i64).unsigned_abs() as u32;
    if len_diff > max_dist {
        return None;
    }

    // Edge cases
    if m == 0 {
        return if n as u32 <= max_dist { Some(n as u32) } else { None };
    }
    if n == 0 {
        return if m as u32 <= max_dist { Some(m as u32) } else { None };
    }

    let max_d = max_dist as usize;

    // Use two rows for space efficiency
    let mut prev = vec![max_dist + 1; n + 1];
    let mut curr = vec![max_dist + 1; n + 1];

    // Initialize first row (within band)
    for (j, val) in prev.iter_mut().enumerate().take(max_d.min(n) + 1) {
        *val = j as u32;
    }

    for i in 1..=m {
        // Band boundaries for this row
        let j_lo = if i > max_d { i - max_d } else { 1 };
        let j_hi = (i + max_d).min(n);

        curr[0] = i as u32;
        let mut row_min = if j_lo == 1 { curr[0] } else { max_dist + 1 };

        for j in j_lo..=j_hi {
            let cost = if a[i - 1] == b[j - 1] { 0u32 } else { 1 };
            let mut val = prev[j - 1] + cost; // substitution
            if prev[j] + 1 < val {
                val = prev[j] + 1; // deletion
            }
            if curr[j - 1] + 1 < val {
                val = curr[j - 1] + 1; // insertion
            }
            curr[j] = val;
            if val < row_min {
                row_min = val;
            }
        }

        // Early termination: entire row exceeds threshold
        if row_min > max_dist {
            return None;
        }

        std::mem::swap(&mut prev, &mut curr);
        // Reset curr for next iteration
        for v in curr.iter_mut() {
            *v = max_dist + 1;
        }
    }

    let dist = prev[n];
    if dist <= max_dist {
        Some(dist)
    } else {
        None
    }
}

fn write_clustered_fasta(path: &str, clusters: &[ClusterResult]) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut w = BufWriter::with_capacity(512 * 1024, file);

    // FASTAptamer format: >Rank-Reads-RPM-Cluster#-ClusterRank-EditDistance
    for cluster in clusters {
        for member in &cluster.members {
            writeln!(
                w,
                ">{}-{}-{:.2}-{}-{}-{}",
                member.original_rank,
                member.count,
                member.rpm,
                cluster.cluster_id,
                member.cluster_rank,
                member.edit_distance,
            )?;
            writeln!(w, "{}", member.sequence)?;
        }
    }

    w.flush()?;
    eprintln!("Clustered FASTA written to: {}", path);
    Ok(())
}

fn write_cluster_summary(path: &str, clusters: &[ClusterResult]) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "cluster_id\tsize\ttotal_reads\tseed_sequence\tseed_count")?;
    for cluster in clusters {
        let seed = &cluster.members[0];
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            cluster.cluster_id,
            cluster.members.len(),
            cluster.total_reads,
            seed.sequence,
            seed.count,
        )?;
    }

    w.flush()?;
    eprintln!("Cluster summary written to: {}", path);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lev_identical() {
        assert_eq!(levenshtein_bounded(b"ATCG", b"ATCG", 3), Some(0));
    }

    #[test]
    fn lev_one_substitution() {
        assert_eq!(levenshtein_bounded(b"ATCG", b"AXCG", 3), Some(1));
    }

    #[test]
    fn lev_one_insertion() {
        assert_eq!(levenshtein_bounded(b"ATCG", b"AXTCG", 3), Some(1));
    }

    #[test]
    fn lev_one_deletion() {
        assert_eq!(levenshtein_bounded(b"ATCG", b"ACG", 3), Some(1));
    }

    #[test]
    fn lev_exceeds_threshold() {
        assert_eq!(levenshtein_bounded(b"AAAA", b"TTTT", 2), None);
    }

    #[test]
    fn lev_at_threshold() {
        assert_eq!(levenshtein_bounded(b"AAAA", b"TTTT", 4), Some(4));
    }

    #[test]
    fn lev_length_diff_exceeds() {
        assert_eq!(levenshtein_bounded(b"AT", b"ATCGATCG", 3), None);
    }

    #[test]
    fn lev_empty() {
        assert_eq!(levenshtein_bounded(b"", b"ATCG", 5), Some(4));
        assert_eq!(levenshtein_bounded(b"ATCG", b"", 5), Some(4));
        assert_eq!(levenshtein_bounded(b"", b"", 0), Some(0));
    }

    #[test]
    fn lev_real_sequences() {
        // 40bp sequences with 3 mismatches
        let a = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        let b = b"ATCGATCGATCXATCGATCXATCGATCGATCXATCGATCG";
        assert_eq!(levenshtein_bounded(a, b, 7), Some(3));
    }

    #[test]
    fn lev_bounded_faster_than_full() {
        // Completely different sequences — should terminate early
        let a = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let b = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        assert_eq!(levenshtein_bounded(a, b, 3), None);
    }
}
