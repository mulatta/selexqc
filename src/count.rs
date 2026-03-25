use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use rayon::prelude::*;

use crate::collector::CountCollector;
use crate::io;

const CHUNK_SIZE: usize = 50_000;

#[derive(Debug, Clone)]
pub enum OutputFormat {
    Parquet,
    Csv,
    Tsv,
    Fasta,
}

pub struct CountConfig {
    pub input: PathBuf,
    pub output_prefix: String,
    pub format: OutputFormat,
    pub rpm: bool,
    pub quiet: bool,
}

pub fn run_count(config: &CountConfig) -> Result<()> {
    // Create output directory if needed
    if let Some(parent) = std::path::Path::new(&config.output_prefix).parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent).ok();
    }

    let mut reader = parse_fastx_file(&config.input).context("Failed to open input file")?;
    let mut count_col = CountCollector::new();

    // Read all sequences into chunks, count in parallel
    let mut chunks: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut current_chunk: Vec<Vec<u8>> = Vec::with_capacity(CHUNK_SIZE);
    let mut total_reads = 0u64;

    while let Some(record) = reader.next() {
        let rec = record.context("Failed to parse record")?;
        current_chunk.push(rec.seq().to_vec());
        total_reads += 1;

        if current_chunk.len() >= CHUNK_SIZE {
            chunks.push(std::mem::take(&mut current_chunk));
            current_chunk = Vec::with_capacity(CHUNK_SIZE);
        }
    }
    if !current_chunk.is_empty() {
        chunks.push(current_chunk);
    }

    // Parallel counting per chunk → merge
    let partial_counts: Vec<ahash::AHashMap<Vec<u8>, u64>> = chunks
        .into_par_iter()
        .map(|chunk| {
            let mut local = ahash::AHashMap::with_capacity(chunk.len() / 2);
            for seq in chunk {
                *local.entry(seq).or_insert(0) += 1;
            }
            local
        })
        .collect();

    // Merge into CountCollector
    for map in partial_counts {
        for (seq, cnt) in map {
            // Use a dummy MatchResult for count-only mode
            count_col.add_count(seq, cnt);
        }
    }
    count_col.set_total_reads(total_reads);

    let entries = count_col.finalize();

    if !config.quiet {
        eprintln!(
            "  {} total reads, {} unique sequences",
            total_reads,
            entries.len()
        );
    }

    // Write output
    match config.format {
        OutputFormat::Parquet => {
            let path = format!("{}.parquet", config.output_prefix);
            let mut meta = HashMap::new();
            meta.insert("selexqc.version".into(), env!("CARGO_PKG_VERSION").into());
            meta.insert("selexqc.input".into(), config.input.display().to_string());
            meta.insert("selexqc.total_reads".into(), total_reads.to_string());
            io::write_count_parquet(&path, &entries, meta, false)?;
        }
        OutputFormat::Csv => {
            let path = format!("{}.csv", config.output_prefix);
            io::write_count_csv(&path, &entries, b',')?;
        }
        OutputFormat::Tsv => {
            let path = format!("{}.tsv", config.output_prefix);
            io::write_count_csv(&path, &entries, b'\t')?;
        }
        OutputFormat::Fasta => {
            let path = format!("{}.fasta", config.output_prefix);
            io::write_ranked_fasta(&path, &entries)?;
        }
    }

    Ok(())
}
