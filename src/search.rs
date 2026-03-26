use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use regex::Regex;

use crate::io;

pub struct SearchConfig {
    pub input: PathBuf,
    pub output_prefix: String,
    pub patterns: Vec<String>,
    pub match_all: bool, // true = AND, false = OR
    pub quiet: bool,
}

pub fn run_search(config: &SearchConfig) -> Result<()> {
    let path = config.input.to_str().unwrap();

    if !config.quiet {
        eprintln!("  Reading: {}", path);
    }

    let entries = io::read_count_parquet(path)?;
    let regexes = compile_iupac_patterns(&config.patterns)?;

    if !config.quiet {
        eprintln!(
            "  Searching {} sequences for {} pattern(s) (mode: {})",
            entries.len(),
            regexes.len(),
            if config.match_all { "AND" } else { "OR" },
        );
    }

    let mut matched = Vec::new();
    for entry in &entries {
        let seq_upper = entry.sequence.to_uppercase();
        let matches: Vec<bool> = regexes.iter().map(|r| r.is_match(&seq_upper)).collect();

        let pass = if config.match_all {
            matches.iter().all(|&m| m)
        } else {
            matches.iter().any(|&m| m)
        };

        if pass {
            matched.push(entry);
        }
    }

    if !config.quiet {
        eprintln!("  {} / {} sequences matched", matched.len(), entries.len());
    }

    // Write matched as ranked FASTA
    let fasta_path = format!("{}.search.fasta", config.output_prefix);
    write_search_fasta(&fasta_path, &matched)?;

    Ok(())
}

fn write_search_fasta(path: &str, entries: &[&io::ParquetEntry]) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut w = BufWriter::new(file);

    for entry in entries {
        writeln!(w, ">{}-{}-{:.2}", entry.rank, entry.count, entry.rpm)?;
        writeln!(w, "{}", entry.sequence)?;
    }

    w.flush()?;
    eprintln!("Search results written to: {}", path);
    Ok(())
}

/// Convert IUPAC ambiguity code patterns to regex.
fn compile_iupac_patterns(patterns: &[String]) -> Result<Vec<Regex>> {
    patterns.iter().map(|p| iupac_to_regex(p)).collect()
}

fn iupac_to_regex(pattern: &str) -> Result<Regex> {
    let mut regex_str = String::with_capacity(pattern.len() * 4);
    for c in pattern.to_uppercase().chars() {
        let expanded = match c {
            'A' => "A",
            'T' | 'U' => "[TU]",
            'G' => "G",
            'C' => "C",
            'R' => "[AG]",
            'Y' => "[CTU]",
            'M' => "[AC]",
            'K' => "[GTU]",
            'S' => "[GC]",
            'W' => "[ATU]",
            'B' => "[CGTU]",
            'D' => "[AGTU]",
            'H' => "[ACTU]",
            'V' => "[ACG]",
            'N' => "[ACGTU]",
            _ => bail!("Invalid IUPAC character '{}' in pattern '{}'", c, pattern),
        };
        regex_str.push_str(expanded);
    }
    Regex::new(&regex_str).with_context(|| format!("Failed to compile pattern '{}'", pattern))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_basic() {
        let re = iupac_to_regex("ATCG").unwrap();
        assert!(re.is_match("ATCG"));
        assert!(!re.is_match("AAAA"));
    }

    #[test]
    fn iupac_ambiguity() {
        let re = iupac_to_regex("RGGARY").unwrap();
        assert!(re.is_match("AGGAGC")); // R=A, A, R=G, Y=C
        assert!(re.is_match("GGGAAT")); // R=G, A, R=A, Y=T
        assert!(!re.is_match("TGGAGC")); // T not in R
    }

    #[test]
    fn iupac_n_matches_all() {
        let re = iupac_to_regex("NNN").unwrap();
        assert!(re.is_match("ATG"));
        assert!(re.is_match("CCC"));
    }

    #[test]
    fn iupac_invalid_char() {
        assert!(iupac_to_regex("ATX").is_err());
    }
}
