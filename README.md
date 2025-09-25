# FastGA-rs

Rust bindings for [FastGA](https://github.com/thegenemyers/FASTGA), a fast genome aligner. No external dependencies - FastGA is compiled and embedded directly.

## Features

- **Zero external dependencies**: FastGA compiled into your Rust binary
- **Process isolation**: Fork/exec architecture prevents memory issues
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

FastGA-rs uses **Rust orchestration with fork/exec isolation**:

1. **Rust API layer** manages configuration and data flow
2. **Fork/exec wrapper** spawns each FastGA utility in separate process
3. **Native FastGA binaries** (FAtoGDB, GIXmake, FastGA) compiled from C

This architecture prevents memory corruption and hanging that occurred with direct FFI.

### Process Flow

```
FASTA files
    ↓
[Fork: FAtoGDB] → Convert to GDB format
    ↓
[Fork: GIXmake] → Build k-mer index
    ↓
[Fork: FastGA]  → Perform alignment
    ↓
PAF output with extended CIGAR
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

## Fork-Based API (Recommended)

For maximum stability, use the fork-based orchestrator:

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