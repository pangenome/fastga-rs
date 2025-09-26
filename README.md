# FastGA-rs

Rust wrapper for [FastGA](https://github.com/thegenemyers/FASTGA), a fast genome aligner. The FastGA binary is compiled and bundled directly - no external dependencies needed.

## Features

- **Zero external dependencies**: FastGA binary bundled with your Rust crate
- **Process isolation**: Subprocess execution prevents memory issues
- **Extended CIGAR**: '=' for matches, 'X' for mismatches
- **Streaming API**: Process large datasets without memory overflow
- **Thread-safe**: Parallel alignment support

## Installation

```toml
[dependencies]
fastga-rs = { git = "https://github.com/pangenome/fastga-rs.git" }
```

## Quick Start

```rust
use fastga_rs::{FastGA, Config};
use std::path::Path;

// Align genomes with default settings (no filtering)
let aligner = FastGA::new(Config::default());
let alignments = aligner.align(
    Path::new("query.fasta"),
    Path::new("target.fasta")
)?;

// Write PAF output
alignments.write_paf("output.paf")?;

// Process alignments
for alignment in alignments.iter() {
    println!("{} -> {}: {:.1}% identity",
        alignment.query_name,
        alignment.target_name,
        alignment.identity() * 100.0);
}
```

## Architecture

FastGA-rs provides a **Rust wrapper around the FastGA binary**:

1. **Rust API layer** - Manages configuration, parameters, and output parsing
2. **Process isolation** - Runs FastGA binary via `std::process::Command`
3. **Native FastGA binary** - Handles everything: FASTA→GDB conversion, indexing, and alignment

This architecture avoids FFI memory issues while providing a clean Rust API.

### Process Flow

```
FASTA files
    ↓
[FastGA binary]
  • Auto-converts to GDB if needed
  • Auto-builds indices if needed
  • Performs alignment
    ↓
PAF output with extended CIGAR
    ↓
[Rust parser]
  • Structured alignment data
  • Iterator APIs
```

## Configuration

```rust
use fastga_rs::Config;

// Default: no filtering
let config = Config::default();

// Custom settings
let config = Config::builder()
    .num_threads(16)
    .min_alignment_length(1000)      // Filter short alignments
    .min_identity(Some(0.90))        // 90% identity threshold
    .output_format(OutputFormat::PafWithX)  // Extended CIGAR
    .keep_intermediates(true)        // Keep .gdb/.gix files
    .build();
```

### Available Parameters

- `num_threads`: Parallelism (default: 8)
- `min_alignment_length`: Length filter (default: 0 - no filter)
- `min_identity`: Identity filter (default: None - no filter)
- `output_format`: PAF variants (default: PafWithX)
- `soft_masking`: Handle lowercase (default: false)
- `keep_intermediates`: Keep temp files (default: false)

## Iteration and Streaming

### Group by Query

Process alignments grouped by query sequence:

```rust
let alignments = aligner.align(query, target)?;
let groups = alignments.group_by_query();

for (query_name, query_alignments) in groups {
    println!("{}: {} alignments", query_name, query_alignments.len());
}
```

### Streaming API

For large-scale alignments, process results without loading everything into memory:

```rust
use fastga_rs::fork_runner::ForkOrchestrator;

let orchestrator = ForkOrchestrator::new(Config::default());

// Each utility runs in isolated process
let paf_output = orchestrator.align(
    Path::new("query.fa"),
    Path::new("target.fa")
)?;
```

## Streaming Large Datasets

Process alignments without loading everything into memory:

```rust
use fastga_rs::simple_runner::run_fastga_streaming;

run_fastga_streaming(
    Path::new("query.fa"),
    Path::new("target.fa"),
    8,      // threads
    0,      // min_length (0 = no filter)
    None,   // min_identity (None = no filter)
    |paf_line| {
        // Process each alignment line
        println!("{}", paf_line);
        true  // Continue (false to stop)
    }
)?;
```

## Building from Source

```bash
git clone https://github.com/pangenome/fastga-rs.git
cd fastga-rs
cargo build --release
cargo test  # Should pass 78/78 tests (100%)
```

## Testing

```bash
# Run all tests
cargo test

# Test with real data
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/chromosomes/chrV.fa.gz
gunzip chrV.fa.gz
cargo run --example align -- chrV.fa chrV.fa > output.paf
# Expect ~1479 alignments for chrV self-alignment
```

## Requirements

- Rust 1.70+
- C compiler (gcc/clang)
- POSIX system (Linux/macOS)

## License

MIT. FastGA included via git subtree under its own license.