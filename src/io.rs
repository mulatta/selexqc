use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{BooleanArray, Float64Array, StringArray, UInt32Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::collector::{
    CombinationCollector, CompositionCollector, ConstructCollector, CountEntry,
};
use crate::construct::Segment;

// --- Parquet I/O ---

pub fn write_count_parquet(
    path: &str,
    entries: &[CountEntry],
    metadata: HashMap<String, String>,
    include_validation: bool,
) -> Result<()> {
    let mut fields = vec![
        Field::new("rank", DataType::UInt64, false),
        Field::new("sequence", DataType::Utf8, false),
        Field::new("count", DataType::UInt64, false),
        Field::new("rpm", DataType::Float64, false),
        Field::new("length", DataType::UInt32, false),
    ];

    if include_validation {
        fields.push(Field::new("is_valid", DataType::Boolean, true));
        fields.push(Field::new("random_region", DataType::Utf8, true));
    }

    let schema = Arc::new(Schema::new(fields));

    let ranks: Vec<u64> = entries.iter().map(|e| e.rank).collect();
    let sequences: Vec<&str> = entries
        .iter()
        .map(|e| std::str::from_utf8(&e.sequence).unwrap_or(""))
        .collect();
    let counts: Vec<u64> = entries.iter().map(|e| e.count).collect();
    let rpms: Vec<f64> = entries.iter().map(|e| e.rpm).collect();
    let lengths: Vec<u32> = entries.iter().map(|e| e.sequence.len() as u32).collect();

    let mut arrays: Vec<Arc<dyn arrow::array::Array>> = vec![
        Arc::new(UInt64Array::from(ranks)),
        Arc::new(StringArray::from(sequences)),
        Arc::new(UInt64Array::from(counts)),
        Arc::new(Float64Array::from(rpms)),
        Arc::new(UInt32Array::from(lengths)),
    ];

    if include_validation {
        // Placeholder: all true for now, will be populated from analyze
        let valids: Vec<bool> = entries.iter().map(|_| true).collect();
        arrays.push(Arc::new(BooleanArray::from(valids)));
        let nulls: Vec<Option<&str>> = entries.iter().map(|_| None).collect();
        arrays.push(Arc::new(StringArray::from(nulls)));
    }

    let batch =
        RecordBatch::try_new(schema.clone(), arrays).context("Failed to create RecordBatch")?;

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;

    let kv_metadata: Vec<parquet::format::KeyValue> = metadata
        .into_iter()
        .map(|(k, v)| parquet::format::KeyValue::new(k, Some(v)))
        .collect();

    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_key_value_metadata(Some(kv_metadata))
        .build();

    let mut writer =
        ArrowWriter::try_new(file, schema, Some(props)).context("Failed to create ArrowWriter")?;
    writer
        .write(&batch)
        .context("Failed to write Parquet data")?;
    writer.close().context("Failed to close Parquet file")?;

    eprintln!("Parquet written to: {}", path);
    Ok(())
}

// --- FASTAptamer-compatible ranked FASTA ---

pub fn write_ranked_fasta(path: &str, entries: &[CountEntry]) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);

    for entry in entries {
        // Format: >RANK-READS-RPM
        write!(writer, ">{}-{}-{:.2}\n", entry.rank, entry.count, entry.rpm)?;
        writer.write_all(&entry.sequence)?;
        writer.write_all(b"\n")?;
    }

    writer.flush()?;
    eprintln!("Ranked FASTA written to: {}", path);
    Ok(())
}

// --- CSV/TSV ---

pub fn write_count_csv(path: &str, entries: &[CountEntry], delimiter: u8) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let writer = BufWriter::with_capacity(512 * 1024, file);
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_writer(writer);

    csv_writer.write_record(["rank", "sequence", "count", "rpm"])?;
    for entry in entries {
        csv_writer.write_record([
            &entry.rank.to_string(),
            std::str::from_utf8(&entry.sequence).unwrap_or(""),
            &entry.count.to_string(),
            &format!("{:.2}", entry.rpm),
        ])?;
    }
    csv_writer.flush()?;
    eprintln!("CSV/TSV written to: {}", path);
    Ok(())
}

// --- Filter Writer ---

pub struct FilterWriter {
    writer: BufWriter<Box<dyn Write + Send>>,
    extract_random: bool,
    random_writer: Option<BufWriter<Box<dyn Write + Send>>>,
}

impl FilterWriter {
    pub fn new(output_prefix: &str, extract_random: bool) -> Result<Self> {
        let path = format!("{}.filtered.fq", output_prefix);
        let file = File::create(&path).with_context(|| format!("Failed to create {}", path))?;
        let writer = BufWriter::with_capacity(1024 * 1024, Box::new(file) as Box<dyn Write + Send>);

        let random_writer = if extract_random {
            let rpath = format!("{}.random.fa", output_prefix);
            let rfile =
                File::create(&rpath).with_context(|| format!("Failed to create {}", rpath))?;
            Some(BufWriter::with_capacity(
                1024 * 1024,
                Box::new(rfile) as Box<dyn Write + Send>,
            ))
        } else {
            None
        };

        Ok(Self {
            writer,
            extract_random,
            random_writer,
        })
    }

    pub fn write_valid(
        &mut self,
        id: &[u8],
        seq: &[u8],
        qual: Option<&[u8]>,
        random_region: Option<&[u8]>,
    ) -> Result<()> {
        // Write as FASTQ if qual available, else FASTA
        if let Some(qual) = qual {
            self.writer.write_all(b"@")?;
            self.writer.write_all(id)?;
            self.writer.write_all(b"\n")?;
            self.writer.write_all(seq)?;
            self.writer.write_all(b"\n+\n")?;
            self.writer.write_all(qual)?;
            self.writer.write_all(b"\n")?;
        } else {
            self.writer.write_all(b">")?;
            self.writer.write_all(id)?;
            self.writer.write_all(b"\n")?;
            self.writer.write_all(seq)?;
            self.writer.write_all(b"\n")?;
        }

        // Write random region
        if self.extract_random {
            if let (Some(writer), Some(region)) = (&mut self.random_writer, random_region) {
                writer.write_all(b">")?;
                writer.write_all(id)?;
                writer.write_all(b"\n")?;
                writer.write_all(region)?;
                writer.write_all(b"\n")?;
            }
        }

        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        if let Some(ref mut w) = self.random_writer {
            w.flush()?;
        }
        Ok(())
    }
}

// --- JSON Report ---

pub fn write_json_report(
    path: &str,
    segments: &[Segment],
    construct_col: &ConstructCollector,
    composition_col: &CompositionCollector,
    combination_col: &CombinationCollector,
    total_reads: u64,
) -> Result<()> {
    let report = build_report_json(
        segments,
        construct_col,
        composition_col,
        combination_col,
        total_reads,
    );

    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    serde_json::to_writer_pretty(file, &report)?;

    eprintln!("JSON report written to: {}", path);
    Ok(())
}

fn build_report_json(
    segments: &[Segment],
    construct_col: &ConstructCollector,
    composition_col: &CompositionCollector,
    combination_col: &CombinationCollector,
    total_reads: u64,
) -> serde_json::Value {
    let mut per_const = Vec::new();
    let mut const_idx = 0;
    for seg in segments {
        if let Segment::Const { sequence } = seg {
            let stats = &construct_col.const_stats()[const_idx];
            per_const.push(serde_json::json!({
                "sequence": String::from_utf8_lossy(sequence),
                "found": stats.found,
                "found_pct": pct(stats.found, total_reads),
                "not_found": stats.not_found,
                "near_miss": {
                    "1mm": stats.near_miss_dist[0],
                    "2mm": stats.near_miss_dist[1],
                    "3mm": stats.near_miss_dist[2],
                    "4mm+": stats.near_miss_dist[3],
                },
            }));
            const_idx += 1;
        }
    }

    let mut per_rand = Vec::new();
    let mut rand_idx = 0;
    for seg in segments {
        if let Segment::Rand { min_len, max_len } = seg {
            let stats = &construct_col.rand_stats()[rand_idx];
            let sorted_dist: std::collections::BTreeMap<usize, u64> =
                stats.length_dist.iter().map(|(&k, &v)| (k, v)).collect();
            per_rand.push(serde_json::json!({
                "expected_range": format!("{}-{}", min_len, max_len),
                "in_range": stats.in_range,
                "in_range_pct": pct(stats.in_range, total_reads),
                "out_of_range": stats.out_of_range,
                "length_distribution": sorted_dist,
            }));
            rand_idx += 1;
        }
    }

    // Composition
    let comp_data: Vec<serde_json::Value> = composition_col
        .per_rand()
        .iter()
        .map(|positions| {
            let pos_data: Vec<serde_json::Value> = positions
                .iter()
                .enumerate()
                .map(|(pos, counts)| {
                    let total: u64 = counts.iter().sum();
                    serde_json::json!({
                        "position": pos + 1,
                        "A": counts[0],
                        "T": counts[1],
                        "G": counts[2],
                        "C": counts[3],
                        "A_pct": pct(counts[0], total),
                        "T_pct": pct(counts[1], total),
                        "G_pct": pct(counts[2], total),
                        "C_pct": pct(counts[3], total),
                    })
                })
                .collect();
            serde_json::json!(pos_data)
        })
        .collect();

    // Combinations (top 20 by count)
    let mut combos: Vec<(u32, u64)> = combination_col
        .combinations()
        .iter()
        .map(|(&mask, &count)| (mask, count))
        .collect();
    combos.sort_by(|a, b| b.1.cmp(&a.1));
    let combo_data: Vec<serde_json::Value> = combos
        .iter()
        .take(20)
        .map(|(mask, count)| {
            let mut pass = Vec::new();
            let mut fail = Vec::new();
            let num_const = combination_col.num_const();
            let num_rand = combination_col.num_rand();
            for i in 0..num_const {
                if mask & (1 << i) != 0 {
                    pass.push(format!("const_{}", i));
                } else {
                    fail.push(format!("const_{}", i));
                }
            }
            for i in 0..num_rand {
                if mask & (1 << (num_const + i)) != 0 {
                    pass.push(format!("rand_{}", i));
                } else {
                    fail.push(format!("rand_{}", i));
                }
            }
            serde_json::json!({
                "pass": pass,
                "fail": fail,
                "count": count,
                "pct": pct(*count, total_reads),
            })
        })
        .collect();

    serde_json::json!({
        "summary": {
            "total_reads": total_reads,
            "valid": construct_col.valid(),
            "valid_pct": pct(construct_col.valid(), total_reads),
            "invalid": total_reads - construct_col.valid(),
        },
        "per_const_segment": per_const,
        "per_rand_segment": per_rand,
        "composition": comp_data,
        "combinations": combo_data,
    })
}

// --- Text Report ---

pub fn write_text_report(
    path: &str,
    segments: &[Segment],
    construct_col: &ConstructCollector,
    combination_col: &CombinationCollector,
    total_reads: u64,
) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "SELEX Construct QC Report")?;
    writeln!(w, "{}", "=".repeat(60))?;
    writeln!(w)?;

    writeln!(w, "Summary:")?;
    writeln!(w, "  Total reads:   {}", total_reads)?;
    writeln!(
        w,
        "  Valid:         {} ({:.2}%)",
        construct_col.valid(),
        pct(construct_col.valid(), total_reads)
    )?;
    writeln!(
        w,
        "  Invalid:       {} ({:.2}%)",
        total_reads - construct_col.valid(),
        pct(total_reads - construct_col.valid(), total_reads)
    )?;
    writeln!(w)?;

    writeln!(w, "Per-Const Segment:")?;
    let mut ci = 0;
    for seg in segments {
        if let Segment::Const { sequence } = seg {
            let stats = &construct_col.const_stats()[ci];
            writeln!(
                w,
                "  [{}] {}: found {} ({:.2}%), near-miss 1mm={} 2mm={} 3mm={}",
                ci,
                String::from_utf8_lossy(sequence),
                stats.found,
                pct(stats.found, total_reads),
                stats.near_miss_dist[0],
                stats.near_miss_dist[1],
                stats.near_miss_dist[2],
            )?;
            ci += 1;
        }
    }
    writeln!(w)?;

    writeln!(w, "Per-Rand Segment:")?;
    let mut ri = 0;
    for seg in segments {
        if let Segment::Rand { min_len, max_len } = seg {
            let stats = &construct_col.rand_stats()[ri];
            writeln!(
                w,
                "  [{}] rand:{}-{}: in_range {} ({:.2}%)",
                ri,
                min_len,
                max_len,
                stats.in_range,
                pct(stats.in_range, total_reads),
            )?;
            ri += 1;
        }
    }
    writeln!(w)?;

    writeln!(w, "Combination Analysis (top 10):")?;
    let mut combos: Vec<(u32, u64)> = combination_col
        .combinations()
        .iter()
        .map(|(&mask, &count)| (mask, count))
        .collect();
    combos.sort_by(|a, b| b.1.cmp(&a.1));
    for (mask, count) in combos.iter().take(10) {
        let all_pass = *mask == combination_col.all_pass_mask();
        let label = if all_pass {
            "ALL PASS".to_string()
        } else {
            format!("mask=0b{:b}", mask)
        };
        writeln!(
            w,
            "  {}: {} ({:.2}%)",
            label,
            count,
            pct(*count, total_reads)
        )?;
    }

    w.flush()?;
    eprintln!("Text report written to: {}", path);
    Ok(())
}

fn pct(n: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        (n as f64 / total as f64) * 100.0
    }
}
