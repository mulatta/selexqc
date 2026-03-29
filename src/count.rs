use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;

use crate::collector::CountCollector;
use crate::io;

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
    if let Some(parent) = std::path::Path::new(&config.output_prefix).parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }

    let (counts, total_reads) = seqtable::count_sequences(&config.input, 0, false)?;
    let count_col = CountCollector::from_counts(counts, total_reads);
    let entries = count_col.finalize();

    if !config.quiet {
        eprintln!(
            "  {} total reads, {} unique sequences",
            total_reads,
            entries.len()
        );
    }

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
