# selexqc

SELEX analysis toolkit: construct QC, sequence counting, and enrichment analysis.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Crates.io](https://img.shields.io/crates/v/selexqc.svg)](https://crates.io/crates/selexqc)

A high-performance replacement for FASTAptamer, written in Rust. Supports RNA/DNA SELEX, Capture-SELEX, and arbitrary library construct definitions.

## Installation

```bash
# From crates.io
cargo install selexqc

# Nix flake
nix profile install github:mulatta/selexqc

# Run without installing
nix run github:mulatta/selexqc -- analyze --help
```

## Quick Start

```bash
# Analyze: validate construct + count sequences (single pass)
selexqc analyze -i round3.fq.gz -o round3 \
  --segment const:GGTAATACGACTCACTATAG \
  --segment rand:35-45 \
  --segment const:ATACCAGCCTCAACTCTGC

# Count only (FASTAptamer-Count replacement)
selexqc count -i round3.fq.gz -o round3 -f parquet --rpm

# Validate only (QC report, no counting)
selexqc validate -i round3.fq.gz -o round3_qc \
  --segment const:ATACCAGCCTCAACTCTGC \
  --segment rand:35-45 \
  --segment const:GGATCCCTCAGCGAGTC
```

## Subcommands

| Command    | Description                                                |
| ---------- | ---------------------------------------------------------- |
| `analyze`  | Validate construct + count in a single pass                |
| `count`    | Count sequence frequencies (FASTAptamer-Count replacement) |
| `validate` | Validate construct structure (no counting)                 |

## Construct Definition

Library constructs are defined with `--segment` flags, in 5' to 3' order:

```bash
--segment const:SEQUENCE   # known fixed region (primer, constant)
--segment rand:MIN-MAX     # variable region (random), length range
```

### Examples

```bash
# Capture-SELEX: 5'primer - random(40bp) - capture constant - 3'primer
--segment const:GGTAATAC --segment rand:35-45 --segment const:ATACCAG --segment const:GGATCC

# Standard SELEX: 5'primer - random - 3'primer
--segment const:GGTAATAC --segment rand:35-45 --segment const:GGATCC

# Simple: just find a constant region
--segment const:ATACCAG
```

## Output Files

### `analyze` outputs

| File                   | Description                                                   |
| ---------------------- | ------------------------------------------------------------- |
| `{prefix}.parquet`     | Sequence count table (rank, sequence, count, rpm)             |
| `{prefix}.fasta`       | FASTAptamer-compatible ranked FASTA (`>RANK-READS-RPM`)       |
| `{prefix}.report.json` | Full QC report (per-segment stats, composition, combinations) |
| `{prefix}.report.txt`  | Human-readable QC summary                                     |
| `{prefix}.filtered.fq` | Valid sequences only (with `--filter`)                        |
| `{prefix}.random.fa`   | Extracted random regions (with `--extract-random`)            |

### `count` outputs

| Format  | Flag         | File               |
| ------- | ------------ | ------------------ |
| Parquet | `-f parquet` | `{prefix}.parquet` |
| CSV     | `-f csv`     | `{prefix}.csv`     |
| TSV     | `-f tsv`     | `{prefix}.tsv`     |
| FASTA   | `-f fasta`   | `{prefix}.fasta`   |

## Key Features

- **Exact match + near-miss reporting**: constant regions matched with SIMD-accelerated exact search (`memchr`). Failed matches analyzed for Hamming distance to report near-misses.
- **Combination analysis**: bitmask-based pass/fail tracking across all segment conditions. Reports UpSet-style intersection statistics.
- **Random region composition**: per-position A/T/G/C counts for each variable region.
- **Parquet-first**: canonical output format with metadata (construct config, version, timestamp). Downstream tools read Parquet directly.
- **FASTAptamer compatible**: ranked FASTA output in `>RANK-READS-RPM` format.
- **Multi-input batch**: process multiple rounds in one command.

## Recommended Pipeline

```bash
# 1. QC/trim with fastp
fastp -i raw.fq.gz -o cleaned.fq.gz

# 2. Analyze with selexqc
selexqc analyze -i cleaned.fq.gz -o results/round3 \
  --segment const:PRIMER5 --segment rand:35-45 --segment const:PRIMER3 \
  --filter --extract-random

# 3. Downstream analysis with Parquet output
# Python: polars.read_parquet("results/round3.parquet")
# SQL:    duckdb "SELECT * FROM 'results/round3.parquet' ORDER BY count DESC LIMIT 20"
```

## Development

### Build from source

```bash
# With Nix
nix develop
cargo build

# Without Nix
cargo build --release
```

### Testing

```bash
cargo test                                        # all tests
cargo run --example generate_fixture              # regenerate test fixtures
```

### Linting

```bash
cargo clippy --all-targets
nix fmt                                           # format nix + rust
```

### Benchmarking

```bash
# Profile with flamegraph (requires nix develop)
cargo flamegraph -- analyze -i large.fq.gz -o bench \
  --segment const:ATCG --segment rand:30-50 --segment const:GCTA -q

# Binary size analysis
cargo bloat --release
```

### Nix

```bash
nix build                   # release build
nix flake check             # run all checks (build + format)
nix run .# -- --version     # run without installing
nix develop                 # enter dev shell
```

## License

MIT
