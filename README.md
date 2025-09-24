# FastGA-rs

A Rust library providing safe, embedded bindings for [FastGA](https://github.com/thegenemyers/FASTGA), a fast genome aligner for high-quality assemblies.

## Features

- **Embedded binaries**: No external dependencies - FastGA is compiled and embedded directly into your Rust binary
- **Extended CIGAR format**: Generates alignments with explicit match ('=') and mismatch ('X') operators
- **Streaming API**: Process alignments as they're generated without storing everything in memory
- **Plane sweep filtering**: Immediate per-query filtering for handling repetitive genomes efficiently
- **Identity scoring**: Calculates alignment identity for downstream filtering algorithms
- **Thread-safe**: Supports parallel alignment operations

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

### Streaming with Immediate Filtering

For repetitive genomes, apply plane sweep filtering immediately to control memory usage:

```rust
use fastga_rs::streaming::align_query_wise_with_sweep;
use fastga_rs::plane_sweep::PlaneSweepConfig;
use fastga_rs::Config;

let sweep_config = PlaneSweepConfig {
    max_per_query: 100,      // Keep top 100 alignments per query
    max_per_target: 1,       // 1:1 mapping per target
    min_identity: 0.90,      // 90% identity threshold
    min_length: 1000,        // Minimum 1kb alignments
    max_overlap: 0.5,        // Filter >50% overlapping alignments
};

let (alignments, stats) = align_query_wise_with_sweep(
    "queries.fa",
    "targets.fa",
    Config::default(),
    sweep_config,
)?;

println!("Kept {} alignments after filtering", alignments.len());
println!("Distribution:\n{}", stats);
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
- **`plane_sweep`**: Filtering algorithms for repetitive sequences
- **`integrated`**: API for embedding in other applications
- **`config`**: Configuration options
- **`alignment`**: Alignment data structures and PAF output

## Requirements

- Rust 1.70 or later
- C compiler (for building FastGA during compilation)
- POSIX-compliant system (Linux, macOS)

## Performance

FastGA-rs provides:
- **Memory efficiency**: Stream processing prevents memory explosion with repetitive genomes
- **Fast alignment**: Leverages FastGA's adaptive seed finding algorithm
- **Parallel processing**: Multi-threaded alignment support
- **Immediate filtering**: Plane sweep filtering during alignment generation

## Use Cases

- **Genome assembly**: Align high-quality assemblies
- **Pangenomics**: All-vs-all alignment of multiple genomes
- **Variant calling**: Extended CIGAR format for precise variant detection
- **Synteny analysis**: Integration with downstream filtering tools

## License

This project is licensed under the MIT License. FastGA is included as a submodule and is subject to its own licensing terms.

## Acknowledgments

- [Gene Myers](https://github.com/thegenemyers) for creating FastGA
- The pangenomics community for feedback and use cases

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## Citation

If you use FastGA-rs in your research, please cite:

```bibtex
@software{fastga-rs,
  author = {Erik Garrison and contributors},
  title = {FastGA-rs: Rust bindings for FastGA genome aligner},
  url = {https://github.com/pangenome/fastga-rs},
  year = {2024}
}
```

And the original FastGA paper (when published).