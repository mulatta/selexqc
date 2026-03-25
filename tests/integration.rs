use std::path::PathBuf;
use tempfile::TempDir;

use selexqc::analyze::{AnalyzeConfig, run_analyze};
use selexqc::construct::Segment;
use selexqc::count::{CountConfig, OutputFormat, run_count};
use selexqc::validate::{ValidateConfig, run_validate};

fn test_segments() -> Vec<Segment> {
    vec![
        Segment::Const {
            sequence: b"AAAAAA".to_vec(),
        },
        Segment::Rand {
            min_len: 8,
            max_len: 12,
        },
        Segment::Const {
            sequence: b"TTTTTT".to_vec(),
        },
    ]
}

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("test.fastq")
}

#[test]
fn analyze_produces_outputs() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("result").to_str().unwrap().to_string();

    run_analyze(&AnalyzeConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        segments: test_segments(),
        strict_boundaries: true,
        near_miss_depth: 3,
        filter: true,
        extract_random: true,
        quiet: true,
    })
    .unwrap();

    // Verify output files exist
    assert!(
        std::path::Path::new(&format!("{}.parquet", output_prefix)).exists(),
        "Parquet file should exist"
    );
    assert!(
        std::path::Path::new(&format!("{}.fasta", output_prefix)).exists(),
        "Ranked FASTA should exist"
    );
    assert!(
        std::path::Path::new(&format!("{}.report.json", output_prefix)).exists(),
        "JSON report should exist"
    );
    assert!(
        std::path::Path::new(&format!("{}.report.txt", output_prefix)).exists(),
        "Text report should exist"
    );
    assert!(
        std::path::Path::new(&format!("{}.filtered.fq", output_prefix)).exists(),
        "Filtered FASTQ should exist"
    );
    assert!(
        std::path::Path::new(&format!("{}.random.fa", output_prefix)).exists(),
        "Random region FASTA should exist"
    );
}

#[test]
fn analyze_json_report_contents() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("result").to_str().unwrap().to_string();

    run_analyze(&AnalyzeConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        segments: test_segments(),
        strict_boundaries: true,
        near_miss_depth: 3,
        filter: false,
        extract_random: false,
        quiet: true,
    })
    .unwrap();

    let json_path = format!("{}.report.json", output_prefix);
    let content = std::fs::read_to_string(&json_path).unwrap();
    let report: serde_json::Value = serde_json::from_str(&content).unwrap();

    // 8 total reads
    assert_eq!(report["summary"]["total_reads"], 8);

    // seq1, seq2, seq3 are valid (perfect match, strict boundaries)
    // seq4: rand too short
    // seq5: missing TTTTTT
    // seq6: missing AAAAAA
    // seq7: leading gap (strict = fail)
    // seq8: duplicate of seq1, valid
    assert_eq!(report["summary"]["valid"], 4);
    assert_eq!(report["summary"]["invalid"], 4);
}

#[test]
fn analyze_ranked_fasta_format() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("result").to_str().unwrap().to_string();

    run_analyze(&AnalyzeConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        segments: test_segments(),
        strict_boundaries: true,
        near_miss_depth: 0,
        filter: false,
        extract_random: false,
        quiet: true,
    })
    .unwrap();

    let fasta = std::fs::read_to_string(format!("{}.fasta", output_prefix)).unwrap();
    let first_line = fasta.lines().next().unwrap();

    // FASTAptamer format: >RANK-READS-RPM
    assert!(
        first_line.starts_with(">1-"),
        "First entry should be rank 1, got: {}",
        first_line
    );
    // Most frequent sequence (seq1 = seq2 = seq8, count=3)
    assert!(
        first_line.contains("-3-"),
        "Most frequent should have count 3, got: {}",
        first_line
    );
}

#[test]
fn count_produces_parquet() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("counts").to_str().unwrap().to_string();

    run_count(&CountConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        format: OutputFormat::Parquet,
        rpm: true,
        quiet: true,
    })
    .unwrap();

    assert!(std::path::Path::new(&format!("{}.parquet", output_prefix)).exists());
}

#[test]
fn count_produces_csv() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("counts").to_str().unwrap().to_string();

    run_count(&CountConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        format: OutputFormat::Csv,
        rpm: true,
        quiet: true,
    })
    .unwrap();

    let csv_path = format!("{}.csv", output_prefix);
    let content = std::fs::read_to_string(&csv_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();

    // Header + 6 unique sequences (8 reads, 3 are duplicates)
    assert_eq!(lines[0], "rank,sequence,count,rpm");
    // seq1=seq2=seq8 (count 3) + seq3 + seq4 + seq5 + seq6 + seq7 = 6 unique
    assert_eq!(lines.len(), 7, "Expected header + 6 unique sequences");
}

#[test]
fn count_produces_fasta() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("counts").to_str().unwrap().to_string();

    run_count(&CountConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        format: OutputFormat::Fasta,
        rpm: false,
        quiet: true,
    })
    .unwrap();

    let fasta_path = format!("{}.fasta", output_prefix);
    let content = std::fs::read_to_string(&fasta_path).unwrap();
    // First header should be the most frequent sequence
    assert!(
        content.starts_with(">1-3-"),
        "Should start with rank 1, count 3"
    );
}

#[test]
fn validate_produces_reports() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("val").to_str().unwrap().to_string();

    run_validate(&ValidateConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        segments: test_segments(),
        strict_boundaries: true,
        near_miss_depth: 3,
        filter: true,
        extract_random: false,
        quiet: true,
    })
    .unwrap();

    assert!(std::path::Path::new(&format!("{}.report.json", output_prefix)).exists());
    assert!(std::path::Path::new(&format!("{}.report.txt", output_prefix)).exists());
    assert!(std::path::Path::new(&format!("{}.filtered.fq", output_prefix)).exists());
}

#[test]
fn validate_filtered_output_only_valid() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("val").to_str().unwrap().to_string();

    run_validate(&ValidateConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        segments: test_segments(),
        strict_boundaries: true,
        near_miss_depth: 0,
        filter: true,
        extract_random: false,
        quiet: true,
    })
    .unwrap();

    let filtered = std::fs::read_to_string(format!("{}.filtered.fq", output_prefix)).unwrap();
    // Each FASTQ record is 4 lines. 4 valid sequences = 16 lines
    let lines: Vec<&str> = filtered.lines().collect();
    assert_eq!(
        lines.len(),
        16,
        "Expected 4 valid sequences (16 lines), got {}",
        lines.len()
    );
}

#[test]
fn analyze_combination_analysis() {
    let tmp = TempDir::new().unwrap();
    let output_prefix = tmp.path().join("result").to_str().unwrap().to_string();

    run_analyze(&AnalyzeConfig {
        input: fixture_path(),
        output_prefix: output_prefix.clone(),
        segments: test_segments(),
        strict_boundaries: true,
        near_miss_depth: 0,
        filter: false,
        extract_random: false,
        quiet: true,
    })
    .unwrap();

    let content = std::fs::read_to_string(format!("{}.report.json", output_prefix)).unwrap();
    let report: serde_json::Value = serde_json::from_str(&content).unwrap();

    // combinations should be an array
    let combos = report["combinations"].as_array().unwrap();
    assert!(!combos.is_empty(), "Should have combination data");

    // First combo should be the most frequent
    let first = &combos[0];
    assert!(first["count"].as_u64().unwrap() > 0);
}
