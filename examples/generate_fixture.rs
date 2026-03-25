//! Generate test FASTQ fixtures for integration tests.
//!
//! Run: cargo run --example generate_fixture

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const CONST_5P: &[u8] = b"AAAAAA";
const CONST_3P: &[u8] = b"TTTTTT";
const QUAL_GOOD: u8 = b'I'; // Q40

fn main() {
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures");
    std::fs::create_dir_all(&out_dir).unwrap();

    let path = out_dir.join("test.fastq");
    let file = File::create(&path).unwrap();
    let mut w = BufWriter::new(file);

    let records = build_test_records();
    for (id, seq, desc) in &records {
        let qual: Vec<u8> = vec![QUAL_GOOD; seq.len()];
        writeln!(w, "@{} {}", id, desc).unwrap();
        w.write_all(seq).unwrap();
        w.write_all(b"\n+\n").unwrap();
        w.write_all(&qual).unwrap();
        w.write_all(b"\n").unwrap();
    }

    w.flush().unwrap();
    eprintln!("Generated {} records → {}", records.len(), path.display());
}

/// Build test records covering all validation scenarios.
///
/// Construct definition: const:AAAAAA - rand:8-12 - const:TTTTTT
fn build_test_records() -> Vec<(&'static str, Vec<u8>, &'static str)> {
    vec![
        // --- Valid sequences ---
        (
            "seq1_valid_perfect",
            build_seq(CONST_5P, b"NNNNNNNNNN", CONST_3P), // rand=10
            "valid:perfect_match",
        ),
        (
            "seq2_valid_duplicate",
            build_seq(CONST_5P, b"NNNNNNNNNN", CONST_3P), // same as seq1
            "valid:duplicate_of_seq1",
        ),
        (
            "seq3_valid_different_rand",
            build_seq(CONST_5P, b"GCATGCATGC", CONST_3P), // rand=10
            "valid:different_random_region",
        ),
        (
            "seq4_valid_min_rand",
            build_seq(CONST_5P, b"NNNNNNNN", CONST_3P), // rand=8 (min)
            "valid:minimum_rand_length",
        ),
        // --- Invalid: rand out of range ---
        (
            "seq5_rand_too_short",
            build_seq(CONST_5P, b"NN", CONST_3P), // rand=2 < 8
            "invalid:rand_too_short",
        ),
        (
            "seq6_rand_too_long",
            build_seq(CONST_5P, b"NNNNNNNNNNNNNN", CONST_3P), // rand=14 > 12
            "invalid:rand_too_long",
        ),
        // --- Invalid: missing const ---
        (
            "seq7_missing_3p",
            [CONST_5P as &[u8], b"NNNNNNNNNN", b"XXXXXX"].concat(),
            "invalid:missing_3p_const",
        ),
        (
            "seq8_missing_5p",
            [b"XXXXXX" as &[u8], b"NNNNNNNNNN", CONST_3P].concat(),
            "invalid:missing_5p_const",
        ),
        // --- Invalid: boundary violations (strict mode) ---
        (
            "seq9_leading_gap",
            [b"GG" as &[u8], CONST_5P, b"NNNNNNNNNN", CONST_3P].concat(),
            "invalid:leading_gap_2bp",
        ),
        (
            "seq10_trailing_gap",
            [CONST_5P as &[u8], b"NNNNNNNNNN", CONST_3P, b"CC"].concat(),
            "invalid:trailing_gap_2bp",
        ),
        // --- Duplicate for count testing ---
        (
            "seq11_third_duplicate",
            build_seq(CONST_5P, b"NNNNNNNNNN", CONST_3P), // same as seq1
            "valid:third_duplicate_of_seq1",
        ),
    ]
}

fn build_seq(prefix: &[u8], middle: &[u8], suffix: &[u8]) -> Vec<u8> {
    let mut seq = Vec::with_capacity(prefix.len() + middle.len() + suffix.len());
    seq.extend_from_slice(prefix);
    seq.extend_from_slice(middle);
    seq.extend_from_slice(suffix);
    seq
}
