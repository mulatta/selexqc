use std::path::PathBuf;

use ahash::AHashMap;
use anyhow::{Context, Result};

use crate::io;

pub struct EnrichConfig {
    pub inputs: Vec<PathBuf>,
    pub output_prefix: String,
    pub labels: Vec<String>,
    pub quiet: bool,
}

pub struct EnrichEntry {
    pub sequence: String,
    pub length: u32,
    /// Per-round: (rank, count, rpm). None if absent in that round.
    pub rounds: Vec<Option<(u64, u64, f64)>>,
    /// Fold enrichment for each adjacent pair: rpm[i+1] / rpm[i].
    pub fold_adjacent: Vec<Option<f64>>,
    /// Fold enrichment first→last: rpm[last] / rpm[first].
    pub fold_overall: Option<f64>,
}

pub fn run_enrich(config: &EnrichConfig) -> Result<()> {
    let paths: Vec<&str> = config.inputs.iter().map(|p| p.to_str().unwrap()).collect();

    // Validate construct consistency
    io::validate_parquet_consistency(&paths)?;

    // Read each round into HashMap<sequence, (rank, count, rpm)>
    let mut round_maps: Vec<AHashMap<String, (u64, u64, f64)>> = Vec::new();
    for (i, path) in paths.iter().enumerate() {
        if !config.quiet {
            eprintln!("  Reading round {}: {}", config.labels[i], path);
        }
        let entries = io::read_count_parquet(path)?;
        let mut map = AHashMap::with_capacity(entries.len());
        for e in entries {
            map.insert(e.sequence, (e.rank, e.count, e.rpm));
        }
        round_maps.push(map);
    }

    // Build union of all sequences
    let mut all_seqs: AHashMap<String, ()> = AHashMap::new();
    for map in &round_maps {
        for seq in map.keys() {
            all_seqs.entry(seq.clone()).or_default();
        }
    }

    // Build enrichment entries
    let num_rounds = round_maps.len();
    let mut entries: Vec<EnrichEntry> = all_seqs
        .keys()
        .map(|seq| {
            let rounds: Vec<Option<(u64, u64, f64)>> = round_maps
                .iter()
                .map(|m| m.get(seq).copied())
                .collect();

            // Adjacent fold enrichment
            let fold_adjacent: Vec<Option<f64>> = (0..num_rounds - 1)
                .map(|i| {
                    match (rounds[i], rounds[i + 1]) {
                        (Some((_, _, rpm_a)), Some((_, _, rpm_b))) if rpm_a > 0.0 => {
                            Some(rpm_b / rpm_a)
                        }
                        _ => None,
                    }
                })
                .collect();

            // Overall fold: last / first
            let fold_overall = match (rounds.first().and_then(|r| *r), rounds.last().and_then(|r| *r)) {
                (Some((_, _, rpm_first)), Some((_, _, rpm_last))) if rpm_first > 0.0 => {
                    Some(rpm_last / rpm_first)
                }
                _ => None,
            };

            EnrichEntry {
                length: seq.len() as u32,
                sequence: seq.clone(),
                rounds,
                fold_adjacent,
                fold_overall,
            }
        })
        .collect();

    // Sort by overall fold enrichment (descending), None last
    entries.sort_by(|a, b| {
        b.fold_overall
            .unwrap_or(f64::NEG_INFINITY)
            .partial_cmp(&a.fold_overall.unwrap_or(f64::NEG_INFINITY))
            .unwrap()
    });

    if !config.quiet {
        eprintln!(
            "  {} unique sequences across {} rounds",
            entries.len(),
            num_rounds
        );
    }

    // Write TSV
    let tsv_path = format!("{}.enrichment.tsv", config.output_prefix);
    write_enrich_tsv(&tsv_path, &entries, &config.labels)?;

    Ok(())
}

fn write_enrich_tsv(
    path: &str,
    entries: &[EnrichEntry],
    labels: &[String],
) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut w = BufWriter::new(file);

    // Header
    write!(w, "sequence\tlength")?;
    for label in labels {
        write!(w, "\t{}_rank\t{}_count\t{}_rpm", label, label, label)?;
    }
    for i in 0..labels.len() - 1 {
        write!(w, "\tfold_{}_vs_{}", labels[i + 1], labels[i])?;
    }
    if labels.len() > 2 {
        write!(w, "\tfold_{}_vs_{}", labels.last().unwrap(), labels.first().unwrap())?;
    }
    writeln!(w)?;

    // Data
    for entry in entries {
        write!(w, "{}\t{}", entry.sequence, entry.length)?;
        for round in &entry.rounds {
            match round {
                Some((rank, count, rpm)) => write!(w, "\t{}\t{}\t{:.2}", rank, count, rpm)?,
                None => write!(w, "\t\t0\t0.00")?,
            }
        }
        for fold in &entry.fold_adjacent {
            match fold {
                Some(f) => write!(w, "\t{:.4}", f)?,
                None => write!(w, "\t")?,
            }
        }
        if labels.len() > 2 {
            match entry.fold_overall {
                Some(f) => write!(w, "\t{:.4}", f)?,
                None => write!(w, "\t")?,
            }
        }
        writeln!(w)?;
    }

    w.flush()?;
    eprintln!("Enrichment TSV written to: {}", path);
    Ok(())
}
