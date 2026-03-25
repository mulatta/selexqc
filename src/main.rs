use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};

use selexqc::construct::{Segment, parse_segment, validate_segments};

#[derive(Parser)]
#[command(name = "selexqc")]
#[command(
    about = "SELEX analysis toolkit: construct QC, sequence counting, and enrichment analysis"
)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Validate construct + count sequences in a single pass
    Analyze {
        /// Input file(s) (FASTQ/FASTA, optionally gzipped)
        #[arg(short, long, required = true, num_args = 1..)]
        input: Vec<PathBuf>,

        /// Output prefix (or directory for multi-input)
        #[arg(short, long)]
        output: String,

        /// Construct segment definition (repeatable, order matters)
        /// Format: "const:SEQUENCE" or "rand:MIN-MAX"
        #[arg(long, required = true, num_args = 1.., value_parser = parse_segment_cli)]
        segment: Vec<Segment>,

        /// Number of threads for parallel processing (0 = auto)
        #[arg(short, long, default_value = "0")]
        threads: usize,

        /// Suppress progress output
        #[arg(short, long)]
        quiet: bool,

        /// Filter: save valid sequences to output file
        #[arg(long)]
        filter: bool,

        /// Extract random regions to separate file
        #[arg(long)]
        extract_random: bool,

        /// Reject sequences with leading/trailing bases outside construct
        #[arg(long, default_value = "true")]
        strict_boundaries: bool,

        /// Near-miss Hamming distance depth for reporting (0 = disable)
        #[arg(long, default_value = "3")]
        near_miss_depth: usize,

        /// Also run enrichment analysis across multiple input rounds
        #[arg(long)]
        enrich: bool,

        /// Include RPM in count output
        #[arg(long)]
        rpm: bool,
    },

    /// Count sequence frequencies (FASTAptamer-Count replacement)
    Count {
        /// Input file(s) (FASTQ/FASTA, optionally gzipped)
        #[arg(short, long, required = true, num_args = 1..)]
        input: Vec<PathBuf>,

        /// Output prefix
        #[arg(short, long)]
        output: String,

        /// Output format
        #[arg(short, long, default_value = "parquet")]
        format: CountFormat,

        /// Include RPM column
        #[arg(long)]
        rpm: bool,

        /// Number of threads (0 = auto)
        #[arg(short, long, default_value = "0")]
        threads: usize,

        /// Suppress progress output
        #[arg(short, long)]
        quiet: bool,
    },

    /// Validate construct structure (no counting)
    Validate {
        /// Input file (FASTQ/FASTA, optionally gzipped)
        #[arg(short, long)]
        input: PathBuf,

        /// Output prefix for report files
        #[arg(short, long)]
        output: String,

        /// Construct segment definition
        #[arg(long, required = true, num_args = 1.., value_parser = parse_segment_cli)]
        segment: Vec<Segment>,

        /// Number of threads (0 = auto)
        #[arg(short, long, default_value = "0")]
        threads: usize,

        /// Suppress progress output
        #[arg(short, long)]
        quiet: bool,

        /// Reject sequences with leading/trailing bases outside construct
        #[arg(long, default_value = "true")]
        strict_boundaries: bool,

        /// Near-miss depth
        #[arg(long, default_value = "3")]
        near_miss_depth: usize,

        /// Filter: save valid sequences
        #[arg(long)]
        filter: bool,

        /// Extract random regions
        #[arg(long)]
        extract_random: bool,
    },
}

#[derive(Debug, Clone, ValueEnum)]
pub enum CountFormat {
    Parquet,
    Csv,
    Tsv,
    Fasta,
}

fn parse_segment_cli(s: &str) -> Result<Segment, String> {
    parse_segment(s).map_err(|e| e.to_string())
}

fn setup_threads(threads: usize) {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok(); // Ignore if already initialized
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Analyze {
            input,
            output,
            segment,
            threads,
            quiet,
            filter,
            extract_random,
            strict_boundaries,
            near_miss_depth,
            enrich: _enrich,
            rpm: _rpm,
        } => {
            validate_segments(&segment)?;
            setup_threads(threads);

            for input_path in &input {
                let out_prefix = if input.len() > 1 {
                    // Multi-input: use output as directory
                    let stem = input_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("output");
                    // Strip .fastq/.fq/.gz suffixes
                    let base = strip_seq_extensions(stem);
                    format!("{}/{}", output, base)
                } else {
                    output.clone()
                };

                if !quiet {
                    eprintln!("Processing: {}", input_path.display());
                }

                selexqc::analyze::run_analyze(&selexqc::analyze::AnalyzeConfig {
                    input: input_path.clone(),
                    output_prefix: out_prefix,
                    segments: segment.clone(),
                    strict_boundaries,
                    near_miss_depth,
                    filter,
                    extract_random,
                    quiet,
                })?;
            }
        }

        Commands::Count {
            input,
            output,
            format,
            rpm,
            threads,
            quiet,
        } => {
            setup_threads(threads);

            for input_path in &input {
                let out_prefix = if input.len() > 1 {
                    let stem = input_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("output");
                    let base = strip_seq_extensions(stem);
                    format!("{}/{}", output, base)
                } else {
                    output.clone()
                };

                if !quiet {
                    eprintln!("Counting: {}", input_path.display());
                }

                selexqc::count::run_count(&selexqc::count::CountConfig {
                    input: input_path.clone(),
                    output_prefix: out_prefix,
                    format: match format {
                        CountFormat::Parquet => selexqc::count::OutputFormat::Parquet,
                        CountFormat::Csv => selexqc::count::OutputFormat::Csv,
                        CountFormat::Tsv => selexqc::count::OutputFormat::Tsv,
                        CountFormat::Fasta => selexqc::count::OutputFormat::Fasta,
                    },
                    rpm,
                    quiet,
                })?;
            }
        }

        Commands::Validate {
            input,
            output,
            segment,
            threads,
            quiet,
            strict_boundaries,
            near_miss_depth,
            filter,
            extract_random,
        } => {
            validate_segments(&segment)?;
            setup_threads(threads);

            if !quiet {
                eprintln!("Validating: {}", input.display());
            }

            selexqc::validate::run_validate(&selexqc::validate::ValidateConfig {
                input,
                output_prefix: output,
                segments: segment,
                strict_boundaries,
                near_miss_depth,
                filter,
                extract_random,
                quiet,
            })?;
        }
    }

    Ok(())
}

fn strip_seq_extensions(name: &str) -> &str {
    let suffixes = [
        ".fastq.gz",
        ".fq.gz",
        ".fasta.gz",
        ".fa.gz",
        ".fastq",
        ".fq",
        ".fasta",
        ".fa",
    ];
    for suffix in &suffixes {
        if let Some(stripped) = name.strip_suffix(suffix) {
            return stripped;
        }
    }
    name
}
