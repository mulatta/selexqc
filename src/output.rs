use anyhow::{Context, Result};
use serde_json::json;
use std::fs::File;
use std::io::Write;

use crate::stats::Stats;
use crate::validate::ValidationMode;

pub fn write_text_report(output_prefix: &str, stats: &Stats) -> Result<()> {
    let filename = format!("{}.validation.txt", output_prefix);
    let mut file = File::create(&filename).context(format!("Failed to create {}", filename))?;

    writeln!(file, "RNA Capture-SELEX Library Validation Report")?;
    writeln!(file, "{}", "=".repeat(70))?;
    writeln!(file)?;

    writeln!(file, "Configuration:")?;
    writeln!(file, "  Constant region: {}", stats.config.constant_region)?;
    writeln!(
        file,
        "  Validation mode: {}",
        if matches!(stats.config.validation_mode, ValidationMode::And) {
            "AND (strict)"
        } else {
            "OR (lenient)"
        }
    )?;

    if let (Some(min), Some(max)) = (stats.config.min_length, stats.config.max_length) {
        writeln!(file, "  Total length range: {} - {} bp", min, max)?;
    }

    if let Some(up_len) = stats.config.upstream_length {
        let tol = stats.config.upstream_tolerance.unwrap_or(0);
        writeln!(
            file,
            "  Expected upstream length: {} bp (+/- {})",
            up_len, tol
        )?;
    }

    if let Some(down_len) = stats.config.downstream_length {
        let tol = stats.config.downstream_tolerance.unwrap_or(0);
        writeln!(
            file,
            "  Expected downstream length: {} bp (+/- {})",
            down_len, tol
        )?;
    }

    if let Some(min_q) = stats.config.min_quality {
        writeln!(file, "  Minimum quality score: {:.1}", min_q)?;
    }
    writeln!(file)?;

    writeln!(file, "Summary Statistics:")?;
    writeln!(file, "  Total sequences: {}", stats.total_sequences)?;
    writeln!(
        file,
        "  Valid sequences: {} ({:.2}%)",
        stats.valid_sequences,
        stats.valid_pct()
    )?;
    writeln!(
        file,
        "  Invalid (filtered) sequences: {} ({:.2}%)",
        stats.invalid_sequences,
        stats.filtered_pct()
    )?;
    writeln!(file)?;

    writeln!(file, "Validation Results:")?;
    writeln!(
        file,
        "  Constant region present: {} ({:.2}%)",
        stats.has_constant_region,
        stats.constant_pct()
    )?;

    if stats.config.min_length.is_some() || stats.config.max_length.is_some() {
        writeln!(
            file,
            "  Correct total length: {} ({:.2}%)",
            stats.correct_total_length,
            stats.correct_length_pct()
        )?;
    }

    if stats.config.upstream_length.is_some() {
        writeln!(
            file,
            "  Correct upstream: {} ({:.2}% of sequences with constant)",
            stats.correct_upstream,
            stats.correct_upstream_pct()
        )?;
    }

    if stats.config.downstream_length.is_some() {
        writeln!(
            file,
            "  Correct downstream: {} ({:.2}% of sequences with constant)",
            stats.correct_downstream,
            stats.correct_downstream_pct()
        )?;
    }

    if stats.config.upstream_length.is_some() && stats.config.downstream_length.is_some() {
        writeln!(
            file,
            "  Correct structure (paired): {} ({:.2}% of sequences with constant)",
            stats.correct_structure_pair,
            stats.correct_structure_pair_pct()
        )?;
    }

    if stats.config.min_quality.is_some() {
        writeln!(
            file,
            "  Correct quality: {} ({:.2}%)",
            stats.correct_quality,
            stats.correct_quality_pct()
        )?;
    }
    writeln!(file)?;

    if !stats.failure_statistics.is_empty() {
        writeln!(
            file,
            "Failure Reasons (sequences can have multiple failures):"
        )?;
        let mut failures: Vec<_> = stats.failure_statistics.iter().collect();
        failures.sort_by_key(|(_, count)| std::cmp::Reverse(*count));

        for (reason, count) in failures {
            let pct = (*count as f64 / stats.invalid_sequences as f64) * 100.0;
            writeln!(
                file,
                "  {}: {} ({:.2}% of invalid sequences)",
                reason, count, pct
            )?;
        }
        writeln!(file)?;
    }

    writeln!(file, "Total Length Distribution:")?;
    write_distribution(&mut file, &stats.length_distribution, stats.total_sequences)?;

    if !stats.upstream_distribution.is_empty() {
        writeln!(file)?;
        writeln!(file, "Upstream Length Distribution (before constant):")?;
        write_distribution(
            &mut file,
            &stats.upstream_distribution,
            stats.has_constant_region,
        )?;
    }

    if !stats.downstream_distribution.is_empty() {
        writeln!(file)?;
        writeln!(file, "Downstream Length Distribution (after constant):")?;
        write_distribution(
            &mut file,
            &stats.downstream_distribution,
            stats.has_constant_region,
        )?;
    }

    if !stats.structure_pair_distribution.is_empty() {
        writeln!(file)?;
        writeln!(file, "Structure Pair Distribution (upstream, downstream):")?;
        write_pair_distribution(
            &mut file,
            &stats.structure_pair_distribution,
            stats.has_constant_region,
        )?;
    }

    eprintln!("Text report written to: {}", filename);
    Ok(())
}

fn write_distribution(
    file: &mut File,
    dist: &std::collections::HashMap<usize, usize>,
    total: usize,
) -> Result<()> {
    let mut items: Vec<_> = dist.iter().collect();
    items.sort_by_key(|(len, _)| *len);

    for (length, count) in items {
        let pct = (*count as f64 / total as f64) * 100.0;
        writeln!(
            file,
            "  {} bp: {:>8} sequences ({:>5.2}%)",
            length, count, pct
        )?;
    }
    Ok(())
}

fn write_pair_distribution(
    file: &mut File,
    dist: &std::collections::HashMap<(usize, usize), usize>,
    total: usize,
) -> Result<()> {
    let mut items: Vec<_> = dist.iter().collect();
    items.sort_by(|(a, _), (b, _)| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    // Show top 20 pairs
    for ((upstream, downstream), count) in items.iter().take(20) {
        let pct = (**count as f64 / total as f64) * 100.0;
        writeln!(
            file,
            "  ({} bp, {} bp): {:>8} sequences ({:>5.2}%)",
            upstream, downstream, count, pct
        )?;
    }

    if items.len() > 20 {
        writeln!(file, "  ... and {} more pairs", items.len() - 20)?;
    }

    Ok(())
}

pub fn write_csv_report(output_prefix: &str, stats: &Stats) -> Result<()> {
    // Total length distribution
    write_csv_dist(
        &format!("{}.length_dist.csv", output_prefix),
        &stats.length_distribution,
        &["length", "count", "percentage"],
        stats.total_sequences,
    )?;

    // Upstream distribution
    if !stats.upstream_distribution.is_empty() {
        write_csv_dist(
            &format!("{}.upstream_dist.csv", output_prefix),
            &stats.upstream_distribution,
            &["upstream_length", "count", "percentage"],
            stats.has_constant_region,
        )?;
    }

    // Downstream distribution
    if !stats.downstream_distribution.is_empty() {
        write_csv_dist(
            &format!("{}.downstream_dist.csv", output_prefix),
            &stats.downstream_distribution,
            &["downstream_length", "count", "percentage"],
            stats.has_constant_region,
        )?;
    }

    // Structure pair distribution
    if !stats.structure_pair_distribution.is_empty() {
        let filename = format!("{}.structure_pairs.csv", output_prefix);
        let mut wtr = csv::Writer::from_path(&filename)?;
        wtr.write_record(&[
            "upstream_length",
            "downstream_length",
            "count",
            "percentage",
        ])?;

        let mut items: Vec<_> = stats.structure_pair_distribution.iter().collect();
        items.sort_by(|(a, _), (b, _)| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        for ((upstream, downstream), count) in items {
            let pct = (*count as f64 / stats.has_constant_region as f64) * 100.0;
            wtr.write_record(&[
                upstream.to_string(),
                downstream.to_string(),
                count.to_string(),
                format!("{:.2}", pct),
            ])?;
        }
        wtr.flush()?;
    }

    eprintln!("CSV reports written to: {}.*.csv", output_prefix);
    Ok(())
}

fn write_csv_dist(
    filename: &str,
    dist: &std::collections::HashMap<usize, usize>,
    headers: &[&str],
    total: usize,
) -> Result<()> {
    let mut wtr = csv::Writer::from_path(filename)?;
    wtr.write_record(headers)?;

    let mut items: Vec<_> = dist.iter().collect();
    items.sort_by_key(|(len, _)| *len);

    for (length, count) in items {
        let pct = (*count as f64 / total as f64) * 100.0;
        wtr.write_record(&[length.to_string(), count.to_string(), format!("{:.2}", pct)])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn write_json_report(output_prefix: &str, stats: &Stats) -> Result<()> {
    let filename = format!("{}.stats.json", output_prefix);
    let file = File::create(&filename).context(format!("Failed to create {}", filename))?;

    serde_json::to_writer_pretty(file, stats)?;
    eprintln!("JSON report written to: {}", filename);
    Ok(())
}

pub fn write_multiqc_report(output_prefix: &str, stats: &Stats) -> Result<()> {
    let filename = format!("{}_mqc.json", output_prefix);

    let validation_mode = if matches!(stats.config.validation_mode, ValidationMode::And) {
        "AND"
    } else {
        "OR"
    };

    let mut general_stats = json!({
        "total_sequences": stats.total_sequences,
        "valid_sequences": stats.valid_sequences,
        "filtered_sequences": stats.invalid_sequences,
        "valid_pct": format!("{:.2}", stats.valid_pct()),
        "filtered_pct": format!("{:.2}", stats.filtered_pct()),
        "has_constant_pct": format!("{:.2}", stats.constant_pct()),
        "validation_mode": validation_mode,
    });

    if stats.config.min_length.is_some() || stats.config.max_length.is_some() {
        general_stats["correct_length_pct"] = json!(format!("{:.2}", stats.correct_length_pct()));
    }

    if stats.config.upstream_length.is_some() && stats.config.downstream_length.is_some() {
        general_stats["correct_structure_pair_pct"] =
            json!(format!("{:.2}", stats.correct_structure_pair_pct()));
    }

    if stats.config.min_quality.is_some() {
        general_stats["correct_quality_pct"] = json!(format!("{:.2}", stats.correct_quality_pct()));
    }

    // Failure reasons
    let mut failure_data = Vec::new();
    for (reason, count) in &stats.failure_statistics {
        failure_data.push(json!({
            "reason": reason.to_string(),
            "count": count,
            "percentage": format!("{:.2}", (*count as f64 / stats.invalid_sequences as f64) * 100.0),
        }));
    }

    let report = json!({
        "id": "slxqc",
        "section_name": "RNA SELEX QC",
        "description": "RNA Capture-SELEX library quality control and filtering statistics",
        "plot_type": "generalstats",
        "pconfig": {
            "namespace": "SELEX QC",
        },
        "data": {
            output_prefix: general_stats
        },
        "failure_reasons": failure_data,
        "filter_efficiency": {
            "total_input": stats.total_sequences,
            "passed_filter": stats.valid_sequences,
            "filtered_out": stats.invalid_sequences,
            "pass_rate": format!("{:.2}", stats.valid_pct()),
        }
    });

    let mut file = File::create(&filename).context(format!("Failed to create {}", filename))?;

    serde_json::to_writer_pretty(&mut file, &report)?;
    eprintln!("MultiQC report written to: {}", filename);
    Ok(())
}
