use std::path::PathBuf;

use anyhow::{Context, Result};
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;

use crate::collector::{CombinationCollector, CompositionCollector, ConstructCollector};
use crate::construct::Segment;
use crate::io;
use crate::matcher::ConstructMatcher;

const CHUNK_SIZE: usize = 10_000;

pub struct ValidateConfig {
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

pub fn run_validate(config: &ValidateConfig) -> Result<()> {
    let matcher = ConstructMatcher::new(&config.segments, config.strict_boundaries);
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

    if let Some(parent) = std::path::Path::new(&config.output_prefix).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).ok();
        }
    }

    let mut reader = parse_fastx_file(&config.input).context("Failed to open input file")?;
    let mut chunk: Vec<RawRecord> = Vec::with_capacity(CHUNK_SIZE);
    let mut total_reads = 0u64;

    while let Some(record) = reader.next() {
        let rec = record.context("Failed to parse record")?;
        chunk.push(RawRecord {
            id: rec.id().to_vec(),
            seq: rec.normalize(false).to_vec(),
            qual: rec.qual().map(|q| q.to_vec()),
        });

        if chunk.len() >= CHUNK_SIZE {
            total_reads += chunk.len() as u64;
            process_chunk(
                &chunk,
                &matcher,
                &mut construct_col,
                &mut comp_col,
                &mut comb_col,
                filter_writer.as_mut(),
            )?;
            chunk.clear();
        }
    }

    if !chunk.is_empty() {
        total_reads += chunk.len() as u64;
        process_chunk(
            &chunk,
            &matcher,
            &mut construct_col,
            &mut comp_col,
            &mut comb_col,
            filter_writer.as_mut(),
        )?;
    }

    if let Some(ref mut fw) = filter_writer {
        fw.flush()?;
    }

    if !config.quiet {
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

    // Write reports
    let json_path = format!("{}.report.json", config.output_prefix);
    let txt_path = format!("{}.report.txt", config.output_prefix);

    io::write_json_report(
        &json_path,
        &config.segments,
        &construct_col,
        &comp_col,
        &comb_col,
        total_reads,
    )?;
    io::write_text_report(
        &txt_path,
        &config.segments,
        &construct_col,
        &comb_col,
        total_reads,
    )?;

    Ok(())
}

fn process_chunk(
    chunk: &[RawRecord],
    matcher: &ConstructMatcher,
    construct_col: &mut ConstructCollector,
    comp_col: &mut CompositionCollector,
    comb_col: &mut CombinationCollector,
    filter_writer: Option<&mut io::FilterWriter>,
) -> Result<()> {
    let match_results: Vec<_> = chunk
        .par_iter()
        .map(|rec| matcher.match_construct(&rec.seq))
        .collect();

    let paired: Vec<(Vec<u8>, crate::matcher::MatchResult)> = chunk
        .iter()
        .zip(match_results.iter())
        .map(|(rec, mr)| (rec.seq.clone(), mr.clone()))
        .collect();

    construct_col.process_chunk(&paired);
    comp_col.process_chunk(&paired);
    comb_col.process_chunk(&paired);

    if let Some(fw) = filter_writer {
        for (rec, result) in chunk.iter().zip(match_results.iter()) {
            if result.matched {
                let random_region = result
                    .rand_regions
                    .first()
                    .and_then(|rm| rm.as_ref().map(|r| &rec.seq[r.start..r.start + r.length]));
                fw.write_valid(&rec.id, &rec.seq, rec.qual.as_deref(), random_region)?;
            }
        }
    }

    Ok(())
}
