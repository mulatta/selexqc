use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;

use crate::collector::{
    CombinationCollector, CompositionCollector, ConstructCollector, CountCollector,
};
use crate::construct::Segment;
use crate::io;
use crate::matcher::ConstructMatcher;

const CHUNK_SIZE: usize = 10_000;

pub struct AnalyzeConfig {
    pub input: PathBuf,
    pub output_prefix: String,
    pub segments: Vec<Segment>,
    pub strict_boundaries: bool,
    pub near_miss_depth: usize,
    pub filter: bool,
    pub extract_random: bool,
    pub quiet: bool,
}

struct RawRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

pub fn run_analyze(config: &AnalyzeConfig) -> Result<()> {
    let matcher = ConstructMatcher::new(&config.segments, config.strict_boundaries);
    let mut count_col = CountCollector::new();
    let mut construct_col = ConstructCollector::new(&config.segments, config.near_miss_depth);
    let mut comp_col = CompositionCollector::new(&config.segments);
    let mut comb_col = CombinationCollector::new(&config.segments);
    let mut filter_writer = if config.filter || config.extract_random {
        Some(io::FilterWriter::new(
            &config.output_prefix,
            config.extract_random,
        )?)
    } else {
        None
    };

    // Create output directory if needed
    if let Some(parent) = std::path::Path::new(&config.output_prefix).parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }

    let mut reader = parse_fastx_file(&config.input).context("Failed to open input file")?;
    let mut chunk: Vec<RawRecord> = Vec::with_capacity(CHUNK_SIZE);

    while let Some(record) = reader.next() {
        let rec = record.context("Failed to parse record")?;
        chunk.push(RawRecord {
            id: rec.id().to_vec(),
            seq: rec.normalize(false).to_vec(),
            qual: rec.qual().map(|q| q.to_vec()),
        });

        if chunk.len() >= CHUNK_SIZE {
            process_chunk(
                &chunk,
                &matcher,
                &mut count_col,
                &mut construct_col,
                &mut comp_col,
                &mut comb_col,
                filter_writer.as_mut(),
            )?;
            chunk.clear();
        }
    }

    // Remaining records
    if !chunk.is_empty() {
        process_chunk(
            &chunk,
            &matcher,
            &mut count_col,
            &mut construct_col,
            &mut comp_col,
            &mut comb_col,
            filter_writer.as_mut(),
        )?;
    }

    // Flush filter writer
    if let Some(ref mut fw) = filter_writer {
        fw.flush()?;
    }
    drop(filter_writer);

    let total_reads = count_col.total_reads();
    let entries = count_col.finalize();

    if !config.quiet {
        eprintln!(
            "  {} total reads, {} unique sequences",
            total_reads,
            entries.len()
        );
        eprintln!(
            "  {}/{} ({:.2}%) valid",
            construct_col.valid(),
            total_reads,
            if total_reads > 0 {
                construct_col.valid() as f64 / total_reads as f64 * 100.0
            } else {
                0.0
            }
        );
    }

    // Build Parquet metadata
    let metadata = build_metadata(config, total_reads);

    // Parallel finalize: write Parquet, ranked FASTA, and reports
    let parquet_path = format!("{}.parquet", config.output_prefix);
    let fasta_path = format!("{}.fasta", config.output_prefix);
    let json_path = format!("{}.report.json", config.output_prefix);
    let txt_path = format!("{}.report.txt", config.output_prefix);

    let errors: Vec<anyhow::Error> = std::thread::scope(|s| {
        let handles: Vec<_> = vec![
            s.spawn(|| io::write_count_parquet(&parquet_path, &entries, metadata, true)),
            s.spawn(|| io::write_ranked_fasta(&fasta_path, &entries)),
            s.spawn(|| {
                io::write_json_report(
                    &json_path,
                    &config.segments,
                    &construct_col,
                    &comp_col,
                    &comb_col,
                    total_reads,
                )
            }),
            s.spawn(|| {
                io::write_text_report(
                    &txt_path,
                    &config.segments,
                    &construct_col,
                    &comb_col,
                    total_reads,
                )
            }),
        ];
        handles
            .into_iter()
            .filter_map(|h| h.join().ok().and_then(|r| r.err()))
            .collect()
    });

    if let Some(first_err) = errors.into_iter().next() {
        return Err(first_err);
    }

    Ok(())
}

fn process_chunk(
    chunk: &[RawRecord],
    matcher: &ConstructMatcher,
    count_col: &mut CountCollector,
    construct_col: &mut ConstructCollector,
    comp_col: &mut CompositionCollector,
    comb_col: &mut CombinationCollector,
    mut filter_writer: Option<&mut io::FilterWriter>,
) -> Result<()> {
    // Parallel validation
    let match_results: Vec<_> = chunk
        .par_iter()
        .map(|rec| matcher.match_construct(&rec.seq))
        .collect();

    // Feed collectors directly (no intermediate clone)
    for (rec, result) in chunk.iter().zip(match_results.iter()) {
        count_col.process_record(&rec.seq);
        construct_col.process_record(&rec.seq, result);
        comp_col.process_record(&rec.seq, result);
        comb_col.process_record(result);

        // Write valid sequences to filter output
        if let Some(fw) = &mut filter_writer
            && result.matched
        {
            let random_region = result
                .rand_regions
                .first()
                .and_then(|rm| rm.as_ref().map(|r| &rec.seq[r.start..r.start + r.length]));
            fw.write_valid(&rec.id, &rec.seq, rec.qual.as_deref(), random_region)?;
        }
    }

    Ok(())
}

fn build_metadata(config: &AnalyzeConfig, total_reads: u64) -> HashMap<String, String> {
    let mut meta = HashMap::new();
    meta.insert("selexqc.version".into(), env!("CARGO_PKG_VERSION").into());
    meta.insert(
        "selexqc.construct".into(),
        serde_json::to_string(&config.segments).unwrap_or_default(),
    );
    meta.insert("selexqc.input".into(), config.input.display().to_string());
    meta.insert(
        "selexqc.timestamp".into(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs().to_string())
            .unwrap_or_default(),
    );
    meta.insert("selexqc.total_reads".into(), total_reads.to_string());
    meta
}
