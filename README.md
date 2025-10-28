# slxqc - RNA Capture-SELEX Library Quality Control

High-performance quality control and filtering tool for RNA Capture-SELEX NGS libraries, written in Rust.

## Features

- ⚡ **Fast**: Parallel processing with zero-copy parsing using `needletail`
- 📁 **Format Support**: FASTA, FASTQ, and gzipped formats (`.fa`, `.fq`, `.fq.gz`)
- 🔍 **Case-Insensitive**: Constant region matching ignores case
- 🏗️ **Structure-Aware**: Validates upstream/downstream regions as paired structures
- 🔀 **Flexible Logic**: AND (strict) or OR (lenient) validation modes
- 🎯 **Filtering**: Save valid sequences to file (FASTA/FASTQ/FASTQ.gz)
- 📊 **MultiQC Integration**: Generate MultiQC-compatible reports
- 🧵 **Scalable**: Multi-threaded processing for large datasets
- 💾 **Memory Efficient**: Streaming processing with minimal memory footprint

## Installation

### From Source (Cargo)

```bash
git clone https://github.com/mulatta/slxqc.git
cd slxqc
cargo build --release
# Binary: target/release/slxqc
```

### Using Nix

```bash
nix build
# Binary: result/bin/slxqc
```

## Quick Start

### Basic Usage

```bash
# Check sequences for constant region presence
slxqc -i library.fq.gz -o results -c TGGCCACCAT
```

### With Structure Validation (Recommended)

```bash
# N40:N10 library (40bp upstream + 10bp constant + 10bp downstream)
slxqc \
  -i library.fq.gz \
  -o results \
  -c TGGCCACCAT \
  --upstream-length 40 \
  --upstream-tolerance 2 \
  --downstream-length 10 \
  --downstream-tolerance 1 \
  --filter \
  --multiqc
```

## Usage

```bash
slxqc [OPTIONS] --input <FILE> --output <PREFIX> --constant <SEQ>
```

### Required Arguments

| Argument                | Description                                           |
| ----------------------- | ----------------------------------------------------- |
| `-i, --input <FILE>`    | Input sequence file (FASTA/FASTQ, optionally gzipped) |
| `-o, --output <PREFIX>` | Output file prefix                                    |
| `-c, --constant <SEQ>`  | Constant region sequence (case-insensitive)           |

### Validation Options

| Option                         | Default | Description                                            |
| ------------------------------ | ------- | ------------------------------------------------------ |
| `--validation-mode <MODE>`     | `and`   | Validation logic: **and** (strict) or **or** (lenient) |
| `--min-length <INT>`           | -       | Minimum total sequence length                          |
| `--max-length <INT>`           | -       | Maximum total sequence length                          |
| `--upstream-length <INT>`      | -       | Expected upstream length (before constant)             |
| `--upstream-tolerance <INT>`   | -       | Upstream length tolerance (+/-)                        |
| `--downstream-length <INT>`    | -       | Expected downstream length (after constant)            |
| `--downstream-tolerance <INT>` | -       | Downstream length tolerance (+/-)                      |
| `-q, --min-quality <FLOAT>`    | -       | Minimum average quality score (FASTQ only)             |

### Filtering Options

| Option                     | Default  | Description                                          |
| -------------------------- | -------- | ---------------------------------------------------- |
| `--filter`                 | disabled | Enable filtering (save valid sequences)              |
| `--filter-format <FORMAT>` | `fasta`  | Output format: **fasta**, **fastq**, or **fastq.gz** |

### Performance & Output

| Option                 | Default    | Description                                    |
| ---------------------- | ---------- | ---------------------------------------------- |
| `-t, --threads <INT>`  | `4`        | Number of threads for parallel processing      |
| `-f, --formats <LIST>` | `txt,json` | Report formats (comma-separated: txt,csv,json) |
| `--multiqc`            | disabled   | Generate MultiQC-compatible report             |

## Validation Logic

### AND Mode (Strict) - Default

**ALL** criteria must pass for a sequence to be valid:

- ✓ Constant region present
- ✓ Total length in range (if specified)
- ✓ Upstream length correct (if specified)
- ✓ Downstream length correct (if specified)
- ✓ Quality sufficient (if specified)

**Use case:** Strict quality control for homogeneous libraries

### OR Mode (Lenient)

**ANY** criterion can pass for a sequence to be valid:

- ✓ Constant region present **OR**
- ✓ Total length in range **OR**
- ✓ Upstream/downstream correct **OR**
- ✓ Quality sufficient

**Use case:** Mixed libraries or exploratory analysis

### Structure Pair Validation

When both `--upstream-length` and `--downstream-length` are specified:

- Validates the **paired structure**: upstream AND downstream together
- Counts sequences where BOTH regions are within tolerance
- Reports paired distribution: `(upstream, downstream)`

**Example:**

```
Expected: 40bp upstream + 10bp constant + 10bp downstream = 60bp total
Tolerance: ±2bp upstream, ±1bp downstream

✓ Valid: 40bp - TGGCCACCAT - 10bp  (exact match)
✓ Valid: 38bp - TGGCCACCAT - 11bp  (within tolerance)
✗ Invalid: 40bp - TGGCCACCAT - 15bp  (downstream too long)
✗ Invalid: 35bp - TGGCCACCAT - 10bp  (upstream too short)
```

## Examples

### Example 1: N40:N10 Library (Strict Validation)

Expected structure: 40bp variable + 10bp constant + 10bp variable = 60bp total

```bash
slxqc \
  -i library.fq.gz \
  -o n40n10_qc \
  -c TGGCCACCAT \
  --validation-mode and \
  --min-length 58 \
  --max-length 62 \
  --upstream-length 40 \
  --upstream-tolerance 2 \
  --downstream-length 10 \
  --downstream-tolerance 1 \
  --filter \
  --filter-format fastq.gz \
  --threads 16 \
  --multiqc
```

**Output files:**

- `n40n10_qc.validation.txt` - Human-readable report
- `n40n10_qc.stats.json` - Complete statistics
- `n40n10_qc.length_dist.csv` - Length distribution
- `n40n10_qc.upstream_dist.csv` - Upstream distribution
- `n40n10_qc.downstream_dist.csv` - Downstream distribution
- `n40n10_qc.structure_pairs.csv` - Paired structure distribution
- `n40n10_qc.filtered.fq.gz` - Valid sequences only
- `n40n10_qc_mqc.json` - MultiQC data

### Example 2: N25:N25 Library

```bash
slxqc \
  -i library.fq.gz \
  -o n25n25_qc \
  -c TGGCCACCAT \
  --upstream-length 25 \
  --upstream-tolerance 2 \
  --downstream-length 25 \
  --downstream-tolerance 2 \
  --filter \
  --multiqc
```

### Example 3: Mixed Library (OR Mode)

Accept sequences with ANY valid characteristic:

```bash
slxqc \
  -i mixed_library.fq.gz \
  -o mixed_qc \
  -c TGGCCACCAT \
  --validation-mode or \
  --min-length 50 \
  --max-length 70
```

### Example 4: Quality Filtering

Filter by quality score (FASTQ only):

```bash
slxqc \
  -i raw.fastq.gz \
  -o qfiltered \
  -c TGGCCACCAT \
  --min-quality 30 \
  --filter \
  --filter-format fastq
```

### Example 5: Structure Analysis Without Filtering

Analyze library structure without saving filtered sequences:

```bash
slxqc \
  -i library.fa \
  -o analysis \
  -c TGGCCACCAT \
  --upstream-length 40 \
  --downstream-length 10

# Review structure pairs
cat analysis.structure_pairs.csv
```

### Example 6: Complete Pipeline with MultiQC

```bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
  slxqc \
    -i ${sample}.fq.gz \
    -o qc/${sample} \
    -c TGGCCACCAT \
    --upstream-length 40 \
    --upstream-tolerance 2 \
    --downstream-length 10 \
    --downstream-tolerance 1 \
    --filter \
    --filter-format fastq.gz \
    --threads 8 \
    --multiqc
done

# Aggregate with MultiQC
multiqc qc/

# Use filtered outputs
nextflow run selex_pipeline \
  --input "qc/*.filtered.fq.gz"
```

## Output Files

### Reports

| File                   | Format | Description                                |
| ---------------------- | ------ | ------------------------------------------ |
| `.validation.txt`      | Text   | Human-readable summary with all statistics |
| `.stats.json`          | JSON   | Complete statistics (machine-readable)     |
| `.length_dist.csv`     | CSV    | Total sequence length distribution         |
| `.upstream_dist.csv`   | CSV    | Upstream region length distribution        |
| `.downstream_dist.csv` | CSV    | Downstream region length distribution      |
| `.structure_pairs.csv` | CSV    | Paired (upstream, downstream) distribution |
| `_mqc.json`            | JSON   | MultiQC-compatible data                    |

### Filtered Sequences

| File              | Format   | Description                                 |
| ----------------- | -------- | ------------------------------------------- |
| `.filtered.fa`    | FASTA    | Valid sequences (uncompressed)              |
| `.filtered.fq`    | FASTQ    | Valid sequences with quality (uncompressed) |
| `.filtered.fq.gz` | FASTQ.gz | Valid sequences with quality (compressed)   |

### MultiQC Report Contents

- Total/valid/filtered sequence counts
- Pass/filter rates
- Validation mode (AND/OR)
- Constant region detection rate
- Structure validation statistics
- Failure reason breakdown
- Filter efficiency metrics

## Understanding Results

### Text Report Example

```
RNA Capture-SELEX Library Validation Report
======================================================================

Configuration:
  Constant region: TGGCCACCAT
  Validation mode: AND (strict)
  Total length range: 58 - 62 bp
  Expected upstream length: 40 bp (+/- 2)
  Expected downstream length: 10 bp (+/- 1)

Summary Statistics:
  Total sequences: 10471867
  Valid sequences: 9525680 (90.96%)
  Invalid (filtered) sequences: 946187 (9.04%)

Validation Results:
  Constant region present: 10450000 (99.79%)
  Correct total length: 10300000 (98.36%)
  Correct upstream: 9800000 (93.78% of sequences with constant)
  Correct downstream: 10200000 (97.61% of sequences with constant)
  Correct structure (paired): 9525680 (91.15% of sequences with constant)

Failure Reasons:
  Incorrect structure: 850320 (89.87% of invalid sequences)
  Incorrect total length: 171867 (18.16% of invalid sequences)
  Missing constant region: 21867 (2.31% of invalid sequences)

Structure Pair Distribution:
  (40, 10): 9000000 sequences (86.12%)
  (39, 11): 300000 sequences ( 2.87%)
  (41,  9): 225680 sequences ( 2.16%)
  ...
```

### Structure Pair Distribution

Shows how many sequences have each (upstream, downstream) combination:

```csv
upstream_length,downstream_length,count,percentage
40,10,9000000,86.12
39,11,300000,2.87
41,9,225680,2.16
38,12,150000,1.43
```

This helps identify:

- Most common library structure
- Distribution of variants
- Potential issues with library preparation

## Performance

<!-- Benchmarks on modern hardware (Intel Xeon, 16 cores): -->

<!-- | Dataset Size | Threads | Processing Time | Memory Usage | Throughput | -->
<!-- | ------------ | ------- | --------------- | ------------ | ---------- | -->
<!-- | 100K reads   | 4       | ~1 second       | ~50 MB       | 100K/s     | -->
<!-- | 1M reads     | 8       | ~8 seconds      | ~200 MB      | 125K/s     | -->
<!-- | 10M reads    | 16      | ~75 seconds     | ~500 MB      | 133K/s     | -->
<!-- | 50M reads    | 32      | ~350 seconds    | ~1 GB        | 143K/s     | -->

**Optimizations:**

- Zero-copy parsing with `needletail`
- Streaming processing (10K sequence chunks)
- Lock-free atomic counters
- Fast substring search with Boyer-Moore (via `memchr`)

## Tips & Best Practices

### Determining Optimal Parameters

1. **Initial exploration** (no structure validation):

```bash
slxqc -i library.fq.gz -o explore -c TGGCCACCAT
cat explore.validation.txt  # Review distributions
```

2. **Check distributions** in the report to see actual upstream/downstream lengths

3. **Set tolerances** based on observed distribution:
   - 95% within range → use that as tolerance
   - Large variation → consider OR mode or wider tolerance

### Recommended Workflows

**For homogeneous libraries (N40:N10, N25:N25):**

- Use AND mode (strict)
- Set structure validation
- Enable filtering
- Use filtered output for downstream analysis

**For mixed or exploratory libraries:**

- Start with OR mode
- Review structure pairs distribution
- Refine parameters based on results

**For quality control in pipelines:**

- Enable MultiQC output
- Use consistent parameters across samples
- Compare QC metrics between samples

### Batch Processing

```bash
# Process all samples in parallel (GNU parallel)
ls *.fq.gz | parallel -j 4 \
  'slxqc -i {} -o qc/{/.} -c TGGCCACCAT \
   --upstream-length 40 --downstream-length 10 \
   --filter --multiqc'

# Aggregate results
multiqc qc/
```

## Troubleshooting

### No Valid Sequences Found

**Check:**

1. Constant region sequence is correct
2. Constant region is in the same format (DNA/RNA, not complement)
3. Review failure statistics in `.validation.txt`

**Debug:**

```bash
# Verify constant region is present
grep -c "TGGCCACCAT" library.fa

# Check if case-sensitive issue (shouldn't be, but verify)
grep -i "tggccaccat" library.fa | head
```

### Low Structure Match Rate

**If paired structure validation fails:**

1. Check `structure_pairs.csv` for actual distribution
2. If distribution shows near-misses, increase tolerances
3. Consider if library actually has mixed structures (use OR mode)

**Example:**

```bash
# Check top structure pairs
head -20 results.structure_pairs.csv

# If you see (38,12) and (42,8) frequently, consider:
# - Increasing tolerance, OR
# - Using OR mode, OR
# - Library has multiple intended structures
```

### Memory Issues

For very large files (>100M reads):

```bash
# Reduce thread count to lower memory usage
slxqc -i huge.fq.gz -o output -c CONST --threads 4

# Or process in chunks (external tool)
split -l 40000000 huge.fq huge_chunk_
for chunk in huge_chunk_*; do
  slxqc -i $chunk -o qc/$(basename $chunk) -c CONST
done
```

## Technical Details

### Dependencies

- `needletail` - Fast FASTA/FASTQ parsing
- `rayon` - Data parallelism
- `memchr` - Fast substring search (Boyer-Moore-Horspool)
- `flate2` - Gzip compression/decompression
- `serde` / `serde_json` - Serialization
- `csv` - CSV writing

### Algorithms

**Constant Region Search:**

- Boyer-Moore-Horspool via `memchr::memmem`
- O(n) average case, O(nm) worst case
- SIMD-accelerated when available

**Parallel Processing:**

- Chunk-based processing (10K sequences per chunk)
- Rayon work-stealing thread pool
- Atomic counters for lock-free statistics
- Minimal memory overhead

**Quality Calculation (FASTQ):**

- Phred quality: Q = -10 \* log10(P)
- ASCII conversion: Q = ASCII - 33
- Average across all positions

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

MIT License - see LICENSE file for details

## Citation

If you use `slxqc` in your research, please cite:

```
slxqc: High-performance quality control for RNA Capture-SELEX libraries
https://github.com/mulatta/slxqc
```

## Support

- **Issues**: https://github.com/mulatta/slxqc/issues
- **Documentation**: https://github.com/mulatta/slxqc/wiki
- **Discussions**: https://github.com/mulatta/slxqc/discussions

## Changelog

### v0.1.0 (Initial Release)

- Fast parallel processing of FASTA/FASTQ/FASTQ.gz
- Case-insensitive constant region detection
- Structure validation (upstream/downstream pairs)
- AND/OR validation modes
- Filtering with multiple output formats
- MultiQC integration
- Comprehensive reporting (TXT/CSV/JSON)
