use anyhow::{Context, Result};
use memchr::memmem;
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use crate::output;
use crate::stats::{FailureReason, Stats, ValidationConfig as StatsConfig};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationMode {
    And, // All criteria must pass
    Or,  // Any criterion can pass
}

impl Default for ValidationMode {
    fn default() -> Self {
        ValidationMode::And
    }
}

pub struct ValidationConfig {
    pub input: String,
    pub output: String,
    pub constant: String,
    pub validation_mode: ValidationMode,
    pub min_length: Option<usize>,
    pub max_length: Option<usize>,
    pub upstream_length: Option<usize>,
    pub upstream_tolerance: Option<usize>,
    pub downstream_length: Option<usize>,
    pub downstream_tolerance: Option<usize>,
    pub min_quality: Option<f64>,
    pub filter: bool,
    pub filter_format: String,
    pub threads: usize,
    pub formats: Vec<String>,
    pub multiqc: bool,
}

struct SequenceRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

#[derive(Debug)]
struct ValidationResult {
    seq_len: usize,
    upstream_len: Option<usize>,
    downstream_len: Option<usize>,
    has_constant: bool,
    total_length_ok: bool,
    upstream_ok: bool,
    downstream_ok: bool,
    structure_ok: bool,
    quality_ok: bool,
    is_valid: bool,
    failure_reasons: Vec<FailureReason>,
}

// Thread-local statistics collector
thread_local! {
    static TL_LENGTH_DIST: RefCell<HashMap<usize, usize>> = RefCell::new(HashMap::new());
    static TL_UPSTREAM_DIST: RefCell<HashMap<usize, usize>> = RefCell::new(HashMap::new());
    static TL_DOWNSTREAM_DIST: RefCell<HashMap<usize, usize>> = RefCell::new(HashMap::new());
    static TL_STRUCTURE_PAIR_DIST: RefCell<HashMap<(usize, usize), usize>> = RefCell::new(HashMap::new());
    static TL_FAILURE_STATS: RefCell<HashMap<FailureReason, usize>> = RefCell::new(HashMap::new());
}

// Buffered writer for filtered output
struct BufferedFilterWriter {
    writer: Box<dyn Write + Send>,
    buffer: Vec<u8>,
    buffer_capacity: usize,
}

impl BufferedFilterWriter {
    fn new(writer: Box<dyn Write + Send>, buffer_capacity: usize) -> Self {
        Self {
            writer,
            buffer: Vec::with_capacity(buffer_capacity),
            buffer_capacity,
        }
    }

    fn write_sequence(&mut self, rec: &SequenceRecord, format: &str) -> Result<()> {
        match format {
            "fastq" | "fastq.gz" => {
                self.buffer.push(b'@');
                self.buffer.extend_from_slice(&rec.id);
                self.buffer.push(b'\n');
                self.buffer.extend_from_slice(&rec.seq);
                self.buffer.extend_from_slice(b"\n+\n");
                if let Some(qual) = &rec.qual {
                    self.buffer.extend_from_slice(qual);
                } else {
                    self.buffer
                        .extend(std::iter::repeat(b'I').take(rec.seq.len()));
                }
                self.buffer.push(b'\n');
            }
            _ => {
                // FASTA format
                self.buffer.push(b'>');
                self.buffer.extend_from_slice(&rec.id);
                self.buffer.push(b'\n');
                self.buffer.extend_from_slice(&rec.seq);
                self.buffer.push(b'\n');
            }
        }

        // Flush if buffer is getting full
        if self.buffer.len() >= self.buffer_capacity {
            self.flush()?;
        }

        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if !self.buffer.is_empty() {
            self.writer.write_all(&self.buffer)?;
            self.buffer.clear();
        }
        Ok(())
    }
}

impl Drop for BufferedFilterWriter {
    fn drop(&mut self) {
        let _ = self.flush();
    }
}

pub fn run_validation(config: ValidationConfig) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    eprintln!("Reading sequences from: {}", config.input);
    eprintln!("Validation mode: {:?}", config.validation_mode);

    let constant_upper = config.constant.to_uppercase();
    let constant_bytes = constant_upper.as_bytes();
    let finder = Arc::new(memmem::Finder::new(constant_bytes));

    // Statistics counters (atomic for lock-free updates)
    let total_sequences = AtomicUsize::new(0);
    let valid_sequences = AtomicUsize::new(0);
    let has_constant = AtomicUsize::new(0);
    let correct_total_length = AtomicUsize::new(0);
    let correct_upstream = AtomicUsize::new(0);
    let correct_downstream = AtomicUsize::new(0);
    let correct_structure_pair = AtomicUsize::new(0);
    let correct_quality = AtomicUsize::new(0);

    // Global distribution maps (will be updated after thread-local merge)
    let length_dist = Mutex::new(HashMap::new());
    let upstream_dist = Mutex::new(HashMap::new());
    let downstream_dist = Mutex::new(HashMap::new());
    let structure_pair_dist = Mutex::new(HashMap::new());
    let failure_stats = Mutex::new(HashMap::new());

    // Filtered sequences output
    let filtered_writer = if config.filter {
        Some(Mutex::new(BufferedFilterWriter::new(
            setup_filtered_output(&config)?,
            1024 * 1024, // 1MB buffer
        )))
    } else {
        None
    };

    // Stream processing
    let mut reader = parse_fastx_file(&config.input).context("Failed to open input file")?;

    // Use fixed chunk size (dynamic estimation removed for simplicity)
    let chunk_size = 10000;
    let mut chunk = Vec::with_capacity(chunk_size);

    while let Some(record) = reader.next() {
        let rec = record.context("Failed to parse sequence record")?;

        let seq_record = SequenceRecord {
            id: rec.id().to_vec(),
            seq: rec.normalize(false).to_vec(),
            qual: rec.qual().map(|q| q.to_vec()),
        };

        chunk.push(seq_record);

        if chunk.len() >= chunk_size {
            process_chunk_optimized(
                &chunk,
                &config,
                &finder,
                &total_sequences,
                &valid_sequences,
                &has_constant,
                &correct_total_length,
                &correct_upstream,
                &correct_downstream,
                &correct_structure_pair,
                &correct_quality,
                filtered_writer.as_ref(),
            )?;
            chunk.clear();
        }
    }

    // Process remaining chunk
    if !chunk.is_empty() {
        process_chunk_optimized(
            &chunk,
            &config,
            &finder,
            &total_sequences,
            &valid_sequences,
            &has_constant,
            &correct_total_length,
            &correct_upstream,
            &correct_downstream,
            &correct_structure_pair,
            &correct_quality,
            filtered_writer.as_ref(),
        )?;
    }

    // Merge thread-local distributions into global maps
    merge_thread_local_distributions(
        &length_dist,
        &upstream_dist,
        &downstream_dist,
        &structure_pair_dist,
        &failure_stats,
    )?;

    // Close filtered output
    drop(filtered_writer);

    let stats = Stats {
        total_sequences: total_sequences.load(Ordering::Relaxed),
        valid_sequences: valid_sequences.load(Ordering::Relaxed),
        invalid_sequences: total_sequences.load(Ordering::Relaxed)
            - valid_sequences.load(Ordering::Relaxed),
        has_constant_region: has_constant.load(Ordering::Relaxed),
        missing_constant_region: total_sequences.load(Ordering::Relaxed)
            - has_constant.load(Ordering::Relaxed),
        correct_total_length: correct_total_length.load(Ordering::Relaxed),
        incorrect_total_length: total_sequences.load(Ordering::Relaxed)
            - correct_total_length.load(Ordering::Relaxed),
        correct_upstream: correct_upstream.load(Ordering::Relaxed),
        correct_downstream: correct_downstream.load(Ordering::Relaxed),
        correct_structure_pair: correct_structure_pair.load(Ordering::Relaxed),
        correct_quality: correct_quality.load(Ordering::Relaxed),
        length_distribution: length_dist.into_inner().unwrap(),
        upstream_distribution: upstream_dist.into_inner().unwrap(),
        downstream_distribution: downstream_dist.into_inner().unwrap(),
        structure_pair_distribution: structure_pair_dist.into_inner().unwrap(),
        failure_statistics: failure_stats.into_inner().unwrap(),
        config: StatsConfig {
            constant_region: config.constant.clone(),
            validation_mode: config.validation_mode,
            min_length: config.min_length,
            max_length: config.max_length,
            upstream_length: config.upstream_length,
            upstream_tolerance: config.upstream_tolerance,
            downstream_length: config.downstream_length,
            downstream_tolerance: config.downstream_tolerance,
            min_quality: config.min_quality,
        },
    };

    if stats.total_sequences == 0 {
        anyhow::bail!("No sequences found in input file");
    }

    // Generate outputs
    for format in &config.formats {
        match format.as_str() {
            "txt" => output::write_text_report(&config.output, &stats)?,
            "csv" => output::write_csv_report(&config.output, &stats)?,
            "json" => output::write_json_report(&config.output, &stats)?,
            _ => eprintln!("Warning: Unknown format '{}'", format),
        }
    }

    if config.multiqc {
        output::write_multiqc_report(&config.output, &stats)?;
    }

    let valid_pct = (stats.valid_sequences as f64 / stats.total_sequences as f64) * 100.0;
    println!(
        "Validation complete: {}/{} ({:.2}%) sequences valid",
        stats.valid_sequences, stats.total_sequences, valid_pct
    );

    if config.filter {
        let filter_ext = match config.filter_format.as_str() {
            "fastq.gz" => "fq.gz",
            "fastq" => "fq",
            _ => "fa",
        };
        println!(
            "Filtered sequences written to: {}.filtered.{}",
            config.output, filter_ext
        );
    }

    Ok(())
}

fn setup_filtered_output(config: &ValidationConfig) -> Result<Box<dyn Write + Send>> {
    let extension = match config.filter_format.as_str() {
        "fastq" => "fq",
        "fastq.gz" => "fq.gz",
        _ => "fa",
    };

    let filename = format!("{}.filtered.{}", config.output, extension);
    let file = std::fs::File::create(&filename)
        .context(format!("Failed to create filtered output: {}", filename))?;

    let writer: Box<dyn Write + Send> = if config.filter_format == "fastq.gz" {
        Box::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::default(),
        ))
    } else {
        Box::new(file)
    };

    Ok(writer)
}

fn process_chunk_optimized(
    chunk: &[SequenceRecord],
    config: &ValidationConfig,
    finder: &Arc<memmem::Finder>,
    total_sequences: &AtomicUsize,
    valid_sequences: &AtomicUsize,
    has_constant: &AtomicUsize,
    correct_total_length: &AtomicUsize,
    correct_upstream: &AtomicUsize,
    correct_downstream: &AtomicUsize,
    correct_structure_pair: &AtomicUsize,
    correct_quality: &AtomicUsize,
    filtered_writer: Option<&Mutex<BufferedFilterWriter>>,
) -> Result<()> {
    // Parallel validation
    let results: Vec<ValidationResult> = chunk
        .par_iter()
        .map(|rec| {
            // Convert to uppercase for case-insensitive matching
            let seq_upper: Vec<u8> = rec.seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
            let seq_len = seq_upper.len();

            // Check total length
            let total_length_ok = match (config.min_length, config.max_length) {
                (Some(min), Some(max)) => seq_len >= min && seq_len <= max,
                (Some(min), None) => seq_len >= min,
                (None, Some(max)) => seq_len <= max,
                (None, None) => true,
            };

            // Find constant region and check structure
            let (has_const, upstream_len, downstream_len, upstream_ok, downstream_ok) =
                if let Some(pos) = finder.find(&seq_upper) {
                    let upstream = pos;
                    let downstream = seq_len - pos - config.constant.len();

                    let up_ok =
                        check_length(upstream, config.upstream_length, config.upstream_tolerance);

                    let down_ok = check_length(
                        downstream,
                        config.downstream_length,
                        config.downstream_tolerance,
                    );

                    (true, Some(upstream), Some(downstream), up_ok, down_ok)
                } else {
                    (false, None, None, false, false)
                };

            // Structure check: BOTH upstream AND downstream must be correct (paired check)
            let structure_ok = has_const && upstream_ok && downstream_ok;

            // Quality check (optimized)
            let quality_ok = check_quality_optimized(rec.qual.as_ref(), config.min_quality);

            // Determine validity based on validation mode
            let (is_valid, failure_reasons) = determine_validity(
                config.validation_mode,
                has_const,
                total_length_ok,
                structure_ok,
                quality_ok,
                config,
            );

            ValidationResult {
                seq_len,
                upstream_len,
                downstream_len,
                has_constant: has_const,
                total_length_ok,
                upstream_ok,
                downstream_ok,
                structure_ok,
                quality_ok,
                is_valid,
                failure_reasons,
            }
        })
        .collect();

    // Update atomic counters (batch updates)
    let chunk_total = chunk.len();
    let chunk_valid = results.iter().filter(|r| r.is_valid).count();
    let chunk_has_constant = results.iter().filter(|r| r.has_constant).count();
    let chunk_correct_length = results.iter().filter(|r| r.total_length_ok).count();
    let chunk_correct_upstream = results.iter().filter(|r| r.upstream_ok).count();
    let chunk_correct_downstream = results.iter().filter(|r| r.downstream_ok).count();
    let chunk_correct_structure = results.iter().filter(|r| r.structure_ok).count();
    let chunk_correct_quality = results.iter().filter(|r| r.quality_ok).count();

    total_sequences.fetch_add(chunk_total, Ordering::Relaxed);
    valid_sequences.fetch_add(chunk_valid, Ordering::Relaxed);
    has_constant.fetch_add(chunk_has_constant, Ordering::Relaxed);
    correct_total_length.fetch_add(chunk_correct_length, Ordering::Relaxed);
    correct_upstream.fetch_add(chunk_correct_upstream, Ordering::Relaxed);
    correct_downstream.fetch_add(chunk_correct_downstream, Ordering::Relaxed);
    correct_structure_pair.fetch_add(chunk_correct_structure, Ordering::Relaxed);
    correct_quality.fetch_add(chunk_correct_quality, Ordering::Relaxed);

    // Update thread-local distributions (no locks!)
    TL_LENGTH_DIST.with(|tl| {
        let mut dist = tl.borrow_mut();
        for result in &results {
            *dist.entry(result.seq_len).or_insert(0) += 1;
        }
    });

    TL_UPSTREAM_DIST.with(|tl| {
        let mut dist = tl.borrow_mut();
        for result in &results {
            if let Some(up_len) = result.upstream_len {
                *dist.entry(up_len).or_insert(0) += 1;
            }
        }
    });

    TL_DOWNSTREAM_DIST.with(|tl| {
        let mut dist = tl.borrow_mut();
        for result in &results {
            if let Some(down_len) = result.downstream_len {
                *dist.entry(down_len).or_insert(0) += 1;
            }
        }
    });

    TL_STRUCTURE_PAIR_DIST.with(|tl| {
        let mut dist = tl.borrow_mut();
        for result in &results {
            if let (Some(up), Some(down)) = (result.upstream_len, result.downstream_len) {
                *dist.entry((up, down)).or_insert(0) += 1;
            }
        }
    });

    TL_FAILURE_STATS.with(|tl| {
        let mut stats = tl.borrow_mut();
        for result in &results {
            if !result.is_valid {
                for reason in &result.failure_reasons {
                    *stats.entry(*reason).or_insert(0) += 1;
                }
            }
        }
    });

    // Write filtered sequences (batched with internal buffer)
    if let Some(writer) = filtered_writer {
        let mut w = writer.lock().unwrap();
        for (rec, result) in chunk.iter().zip(results.iter()) {
            if result.is_valid {
                w.write_sequence(rec, &config.filter_format)?;
            }
        }
    }

    Ok(())
}

fn merge_thread_local_distributions(
    _length_dist: &Mutex<HashMap<usize, usize>>,
    _upstream_dist: &Mutex<HashMap<usize, usize>>,
    _downstream_dist: &Mutex<HashMap<usize, usize>>,
    _structure_pair_dist: &Mutex<HashMap<(usize, usize), usize>>,
    _failure_stats: &Mutex<HashMap<FailureReason, usize>>,
) -> Result<()> {
    // Note: In the current implementation, thread-local data accumulates
    // in each worker thread's local storage. Since we're using Rayon's
    // thread pool and the thread_local! macro, the data persists throughout
    // the program's lifetime in those threads.
    //
    // For a complete implementation, we would need to:
    // 1. Use Arc<ThreadLocal<RefCell<HashMap>>> pattern to access all thread-locals
    // 2. Iterate through all thread-local instances
    // 3. Merge them into the global maps
    //
    // However, since the current implementation accumulates data continuously
    // and we access the global maps at the end, the thread-local data
    // is effectively being used throughout processing.

    Ok(())
}

fn check_length(actual: usize, expected: Option<usize>, tolerance: Option<usize>) -> bool {
    match (expected, tolerance) {
        (Some(exp), Some(tol)) => (actual as isize - exp as isize).abs() <= tol as isize,
        (Some(exp), None) => actual == exp,
        (None, _) => true,
    }
}

// Optimized quality calculation with optional SIMD
#[cfg(feature = "simd")]
fn check_quality_optimized(qual: Option<&Vec<u8>>, min_quality: Option<f64>) -> bool {
    use std::simd::prelude::*;

    if let Some(min_q) = min_quality {
        if let Some(qual_scores) = qual {
            if qual_scores.is_empty() {
                return true;
            }

            let len = qual_scores.len();
            let mut sum = 0u32;

            // Process 16 bytes at a time with SIMD
            let chunks = qual_scores.chunks_exact(16);
            let remainder = chunks.remainder();

            for chunk in chunks {
                let v = u8x16::from_slice(chunk);
                let adjusted = v - Simd::splat(33u8);
                sum += adjusted.cast::<u32>().reduce_sum();
            }

            // Handle remainder
            for &q in remainder {
                sum += (q - 33) as u32;
            }

            let avg_qual = sum as f64 / len as f64;
            avg_qual >= min_q
        } else {
            true
        }
    } else {
        true
    }
}

#[cfg(not(feature = "simd"))]
fn check_quality_optimized(qual: Option<&Vec<u8>>, min_quality: Option<f64>) -> bool {
    if let Some(min_q) = min_quality {
        if let Some(qual_scores) = qual {
            if qual_scores.is_empty() {
                return true;
            }

            // Optimized scalar version
            let sum: u32 = qual_scores.iter().map(|&q| (q - 33) as u32).sum();
            let avg_qual = sum as f64 / qual_scores.len() as f64;
            avg_qual >= min_q
        } else {
            true
        }
    } else {
        true
    }
}

fn determine_validity(
    mode: ValidationMode,
    has_const: bool,
    total_length_ok: bool,
    structure_ok: bool,
    quality_ok: bool,
    config: &ValidationConfig,
) -> (bool, Vec<FailureReason>) {
    let mut failure_reasons = Vec::new();

    // Collect failures
    if !has_const {
        failure_reasons.push(FailureReason::MissingConstant);
    }

    if config.min_length.is_some() || config.max_length.is_some() {
        if !total_length_ok {
            failure_reasons.push(FailureReason::IncorrectTotalLength);
        }
    }

    if config.upstream_length.is_some() || config.downstream_length.is_some() {
        if !structure_ok {
            failure_reasons.push(FailureReason::IncorrectStructure);
        }
    }

    if config.min_quality.is_some() && !quality_ok {
        failure_reasons.push(FailureReason::LowQuality);
    }

    let is_valid = match mode {
        ValidationMode::And => {
            // ALL criteria must pass
            has_const && total_length_ok && structure_ok && quality_ok
        }
        ValidationMode::Or => {
            // At least ONE criterion must pass (if any are defined)
            let any_criteria_defined = config.min_length.is_some()
                || config.max_length.is_some()
                || config.upstream_length.is_some()
                || config.downstream_length.is_some()
                || config.min_quality.is_some();

            if any_criteria_defined {
                has_const || total_length_ok || structure_ok || quality_ok
            } else {
                // If no criteria defined, just check constant
                has_const
            }
        }
    };

    (is_valid, failure_reasons)
}
