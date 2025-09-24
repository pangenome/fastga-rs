# FastGA-rs

A Rust library providing safe, embedded bindings for [FastGA](https://github.com/thegenemyers/FASTGA), a fast genome aligner for high-quality assemblies.

## Features

- **Embedded binaries**: No external dependencies - FastGA is compiled and embedded directly into your Rust binary
- **Extended CIGAR format**: Generates alignments with explicit match ('=') and mismatch ('X') operators
- **Query-complete alignment sets**: Get ALL alignments for each query before processing the next
- **Automatic backpressure**: Prevents memory overflow when processing large datasets
- **Identity scoring**: Calculates alignment identity for downstream filtering algorithms
- **Thread-safe**: Supports parallel alignment operations
- **Memory efficient**: Bounded memory usage with streaming

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

### Query-Complete Alignment Sets

**IMPORTANT**: FastGA processes queries sequentially and completes ALL alignments for each query before moving to the next. This allows you to get complete alignment sets for each query:

```rust
use fastga_rs::{QueryAlignmentIterator, Config};

// Create iterator with backpressure (won't overflow memory)
let query_sets = QueryAlignmentIterator::new(
    "queries.fa",      // Multiple query sequences
    "targets.fa",      // Target database
    Config::default(),
    1,                 // Buffer size (1 = strict backpressure)
)?;

// Process each query's complete alignment set
for query_set in query_sets {
    let query_set = query_set?;

    println!("Query: {} has {} alignments",
             query_set.query_name,
             query_set.alignment_count());

    // All alignments for this query are available
    for alignment in &query_set.alignments {
        println!("  -> {} (identity: {:.1}%)",
                 alignment.target_name,
                 alignment.identity() * 100.0);
    }

    // Get best hit
    if let Some(best) = query_set.best_by_identity() {
        println!("Best hit: {} at {:.1}% identity",
                 best.target_name,
                 best.identity() * 100.0);
    }
}
```

### Simple Query Processing

For simpler use cases, use the high-level API:

```rust
use fastga_rs::{align_queries, Config};

align_queries(
    "queries.fa",
    "targets.fa",
    Config::default(),
    |query_set| {
        // Process complete alignment set for this query
        println!("Processing {} with {} hits",
                 query_set.query_name,
                 query_set.alignment_count());

        // Return true to continue, false to stop
        Ok(true)
    }
)?;
```

### Memory-Efficient Streaming

The query iterator uses bounded channels with backpressure. When you process slowly, FastGA automatically pauses, preventing memory overflow:

```rust
let query_sets = QueryAlignmentIterator::new(
    "huge_queries.fa",
    "huge_database.fa",
    Config::default(),
    1,  // Only buffer 1 query ahead (strict backpressure)
)?;

for query_set in query_sets {
    let query_set = query_set?;

    // FastGA is paused here if we're slow!
    expensive_computation(&query_set)?;

    // FastGA resumes when we loop back
}
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