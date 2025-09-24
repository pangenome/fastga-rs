# FastGA-rs

A Rust library providing safe, embedded bindings for [FastGA](https://github.com/thegenemyers/FASTGA), a fast genome aligner for high-quality assemblies.

## Features

- **Embedded binaries**: No external dependencies - FastGA is compiled and embedded directly into your Rust binary
- **Extended CIGAR format**: Generates alignments with explicit match ('=') and mismatch ('X') operators
- **Streaming API**: Process alignments as they're generated without storing everything in memory
- **Identity scoring**: Calculates alignment identity for downstream filtering algorithms
- **Thread-safe**: Supports parallel alignment operations
- **Memory efficient**: Suitable for large genome alignments

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
fastga-rs = { git = "https://github.com/pangenome/fastga-rs.git" }
```

## Quick Start

```rust
use fastga_rs::{FastGA, Config};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create aligner with default configuration
    let aligner = FastGA::new(Config::default())?;

    // Align two genome files
    let alignments = aligner.align_files(
        Path::new("genome1.fasta"),
        Path::new("genome2.fasta")
    )?;

    // Output in PAF format with extended CIGAR
    println!("{}", alignments.to_paf()?);

    // Access identity scores for filtering
    for alignment in alignments.alignments {
        let identity = alignment.identity();
        let length = alignment.query_end - alignment.query_start;
        let score = identity * (length as f64).ln(); // Plane sweep score
        println!("Alignment: identity={:.2}%, score={:.2}",
                 identity * 100.0, score);
    }

    Ok(())
}
```

## Advanced Usage

### Streaming API

Process alignments as they're generated without loading everything into memory:

```rust
use fastga_rs::streaming::StreamingAligner;
use fastga_rs::Config;

let mut aligner = StreamingAligner::new(Config::default());

// Add filters
aligner
    .filter_min_identity(0.9)
    .filter_min_length(1000);

// Process alignments with a callback
let stats = aligner.align_files(
    "genome1.fa",
    "genome2.fa",
    |alignment| {
        // Process each alignment as it's generated
        println!("Found: {} -> {}",
                 alignment.query_name,
                 alignment.target_name);
        true  // Keep alignment
    }
)?;

println!("Processed {} alignments", stats.total_alignments);
```

### Configuration Presets

```rust
use fastga_rs::{FastGA, Config};

// High sensitivity for divergent sequences
let aligner = FastGA::new(Config::high_sensitivity())?;

// Fast mode for similar sequences
let aligner = FastGA::new(Config::fast())?;

// Optimized for repetitive genomes
let aligner = FastGA::new(Config::repetitive_genomes())?;

// Custom configuration
let config = Config::builder()
    .min_identity(0.85)
    .min_alignment_length(5000)
    .num_threads(8)
    .build();
let aligner = FastGA::new(config)?;
```

### Integration with Other Tools

FastGA-rs is designed to integrate seamlessly with genome analysis pipelines. See [fastga-rs-integration.md](fastga-rs-integration.md) for detailed integration guide with [SweepGA](https://github.com/ekg/sweepga).

## Building from Source

```bash
git clone https://github.com/pangenome/fastga-rs.git
cd fastga-rs

# Initialize submodules (FastGA source)
git submodule update --init --recursive

# Build and run tests
cargo build --release
cargo test
```

## Architecture

The library consists of several modules:

- **`embedded`**: Manages embedded FastGA binaries
- **`streaming`**: Streaming alignment processing
- **`config`**: Configuration options
- **`alignment`**: Alignment data structures and PAF output
- **`error`**: Error handling types

## Requirements

- Rust 1.70 or later
- C compiler (for building FastGA during compilation)
- POSIX-compliant system (Linux, macOS)

## License

This project is licensed under the MIT License. FastGA source is included via git subtree and is subject to its own licensing terms.

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

### Updating FastGA

FastGA is included as a git subtree (not submodule) for crate compatibility. To update to the latest FastGA version:

```bash
git subtree pull --prefix=deps/fastga https://github.com/thegenemyers/FASTGA.git main --squash
```

## Citation

If you use FastGA-rs in your research, please cite the original FastGA:

```bibtex
@article{myers2025fastga,
  author = {Gene Myers and Richard Durbin and Chenxi Zhou},
  title = {FastGA: Fast Genome Alignment},
  journal = {bioRxiv},
  year = {2025},
  doi = {10.1101/2025.06.15.659750}
}
```