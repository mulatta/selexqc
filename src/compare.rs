use std::path::PathBuf;

use ahash::AHashMap;
use anyhow::{Context, Result};

use crate::io;

pub struct CompareConfig {
    pub input_x: PathBuf,
    pub input_y: PathBuf,
    pub output_prefix: String,
    pub label_x: String,
    pub label_y: String,
    pub quiet: bool,
}

pub fn run_compare(config: &CompareConfig) -> Result<()> {
    let path_x = config.input_x.to_str().unwrap();
    let path_y = config.input_y.to_str().unwrap();

    io::validate_parquet_consistency(&[path_x, path_y])?;

    if !config.quiet {
        eprintln!("  Reading {}: {}", config.label_x, path_x);
        eprintln!("  Reading {}: {}", config.label_y, path_y);
    }

    let entries_x = io::read_count_parquet(path_x)?;
    let entries_y = io::read_count_parquet(path_y)?;

    let map_x: AHashMap<String, (u64, f64)> = entries_x
        .into_iter()
        .map(|e| (e.sequence, (e.count, e.rpm)))
        .collect();
    let map_y: AHashMap<String, (u64, f64)> = entries_y
        .into_iter()
        .map(|e| (e.sequence, (e.count, e.rpm)))
        .collect();

    // Classify sequences
    let mut both = Vec::new();
    let mut x_only = 0u64;
    let mut y_only = 0u64;

    for (seq, (count_x, rpm_x)) in &map_x {
        if let Some((count_y, rpm_y)) = map_y.get(seq) {
            let log2_ratio = if *rpm_x > 0.0 {
                (rpm_y / rpm_x).log2()
            } else {
                f64::INFINITY
            };
            both.push(CompareEntry {
                sequence: seq.clone(),
                count_x: *count_x,
                rpm_x: *rpm_x,
                count_y: *count_y,
                rpm_y: *rpm_y,
                log2_ratio,
            });
        } else {
            x_only += 1;
        }
    }
    for seq in map_y.keys() {
        if !map_x.contains_key(seq) {
            y_only += 1;
        }
    }

    // Sort by log2 ratio descending
    both.sort_by(|a, b| {
        b.log2_ratio
            .partial_cmp(&a.log2_ratio)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Histogram: 102 bins (log2 -5 to +5, 0.1 step, plus 2 outlier bins)
    let histogram = build_histogram(&both);

    if !config.quiet {
        eprintln!(
            "  Common: {}, {}-only: {}, {}-only: {}",
            both.len(),
            config.label_x,
            x_only,
            config.label_y,
            y_only
        );
    }

    // Write TSV
    let tsv_path = format!("{}.compare.tsv", config.output_prefix);
    write_compare_tsv(&tsv_path, &both, &config.label_x, &config.label_y)?;

    // Write histogram
    let hist_path = format!("{}.compare.histogram.tsv", config.output_prefix);
    write_histogram_tsv(&hist_path, &histogram)?;

    Ok(())
}

struct CompareEntry {
    sequence: String,
    count_x: u64,
    rpm_x: f64,
    count_y: u64,
    rpm_y: f64,
    log2_ratio: f64,
}

fn build_histogram(entries: &[CompareEntry]) -> Vec<(String, u64)> {
    // 102 bins: [-inf, -5), [-5.0, -4.9), ..., [4.9, 5.0), [5.0, +inf)
    let mut bins = vec![0u64; 102];

    for entry in entries {
        let r = entry.log2_ratio;
        let idx = if r < -5.0 {
            0
        } else if r >= 5.0 {
            101
        } else {
            // Map [-5.0, 5.0) to bins 1..101
            ((r + 5.0) / 0.1).floor() as usize + 1
        };
        bins[idx.min(101)] += 1;
    }

    let mut result = Vec::with_capacity(102);
    result.push(("< -5.0".to_string(), bins[0]));
    for (i, &count) in bins[1..101].iter().enumerate() {
        let low = -5.0 + i as f64 * 0.1;
        result.push((format!("{:.1}", low), count));
    }
    result.push((">= 5.0".to_string(), bins[101]));
    result
}

fn write_compare_tsv(
    path: &str,
    entries: &[CompareEntry],
    label_x: &str,
    label_y: &str,
) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "sequence\t{}_count\t{}_rpm\t{}_count\t{}_rpm\tlog2_ratio",
        label_x, label_x, label_y, label_y
    )?;

    for e in entries {
        writeln!(
            w,
            "{}\t{}\t{:.2}\t{}\t{:.2}\t{:.4}",
            e.sequence, e.count_x, e.rpm_x, e.count_y, e.rpm_y, e.log2_ratio
        )?;
    }

    w.flush()?;
    eprintln!("Compare TSV written to: {}", path);
    Ok(())
}

fn write_histogram_tsv(path: &str, histogram: &[(String, u64)]) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "log2_bin\tcount")?;
    for (label, count) in histogram {
        writeln!(w, "{}\t{}", label, count)?;
    }

    w.flush()?;
    eprintln!("Histogram TSV written to: {}", path);
    Ok(())
}
