use anyhow::Result;
use clap::Parser;

mod output;
mod stats;
mod validate;

#[derive(Parser)]
#[command(name = "selexqc")]
#[command(about = "RNA Capture-SELEX Library Quality Control", long_about = None)]
#[command(version)]
struct Cli {
    /// Input file (FASTA/FASTQ, optionally gzipped)
    #[arg(short, long)]
    input: String,

    /// Output prefix for result files
    #[arg(short, long)]
    output: String,

    /// Constant region sequence to search for
    #[arg(short, long)]
    constant: String,

    /// Validation mode: how to combine multiple criteria
    /// - 'and': All criteria must pass (strict)
    /// - 'or': Any criterion can pass (lenient)
    #[arg(long, default_value = "and")]
    validation_mode: String,

    /// Minimum total sequence length
    #[arg(long)]
    min_length: Option<usize>,

    /// Maximum total sequence length
    #[arg(long)]
    max_length: Option<usize>,

    /// Expected upstream length (before constant region)
    #[arg(long)]
    upstream_length: Option<usize>,

    /// Upstream length tolerance (+/-)
    #[arg(long)]
    upstream_tolerance: Option<usize>,

    /// Expected downstream length (after constant region)
    #[arg(long)]
    downstream_length: Option<usize>,

    /// Downstream length tolerance (+/-)
    #[arg(long)]
    downstream_tolerance: Option<usize>,

    /// Minimum quality score threshold (for FASTQ)
    #[arg(short, long)]
    min_quality: Option<f64>,

    /// Filter and save valid sequences to output file
    #[arg(long)]
    filter: bool,

    /// Output format for filtered sequences (fasta, fastq, fastq.gz)
    #[arg(long, default_value = "fasta")]
    filter_format: String,

    /// Number of threads for parallel processing
    #[arg(short, long, default_value = "4")]
    threads: usize,

    /// Output report formats (txt,csv,json - comma separated)
    #[arg(short = 'f', long, default_value = "txt,json")]
    formats: String,

    /// Generate MultiQC compatible output
    #[arg(long)]
    multiqc: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let validation_mode = match cli.validation_mode.to_lowercase().as_str() {
        "and" => validate::ValidationMode::And,
        "or" => validate::ValidationMode::Or,
        _ => anyhow::bail!("Invalid validation mode. Use 'and' or 'or'"),
    };

    validate::run_validation(validate::ValidationConfig {
        input: cli.input,
        output: cli.output,
        constant: cli.constant,
        validation_mode,
        min_length: cli.min_length,
        max_length: cli.max_length,
        upstream_length: cli.upstream_length,
        upstream_tolerance: cli.upstream_tolerance,
        downstream_length: cli.downstream_length,
        downstream_tolerance: cli.downstream_tolerance,
        min_quality: cli.min_quality,
        filter: cli.filter,
        filter_format: cli.filter_format,
        threads: cli.threads,
        formats: cli
            .formats
            .split(',')
            .map(|s| s.trim().to_string())
            .collect(),
        multiqc: cli.multiqc,
    })?;

    Ok(())
}
