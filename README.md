# FastGA-RS

Rust bindings for [FastGA](https://github.com/thegenemyers/FASTGA), a fast genome aligner for high-quality genome assemblies.

## Features

- **High Performance**: Leverages FastGA's novel adaptive seed finding algorithm
- **Extended CIGAR Support**: Produces detailed alignments with explicit match ('=') and mismatch ('X') operators
- **Flexible API**: Both simple file-based and advanced streaming interfaces
- **Memory Efficient**: Uses trace point compression for alignment storage
- **Thread-Safe**: Supports parallel alignment operations
- **Configurable**: Fine-tune alignment parameters for different use cases

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
fastga-rs = "0.1"
```

### System Requirements

FastGA-RS requires the following system libraries:
- pthread
- zlib
- A C compiler (gcc or clang)

On Ubuntu/Debian:
```bash
sudo apt-get install build-essential zlib1g-dev
```

On macOS:
```bash
# Install Xcode command line tools if not already installed
xcode-select --install
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

    // Process alignments
    for alignment in alignments.iter() {
        println!("{} -> {}: {:.1}% identity",
            alignment.query_name,
            alignment.target_name,
            alignment.identity() * 100.0
        );
    }

    // Output in PAF format with extended CIGAR
    let paf_output = alignments.to_paf()?;
    println!("{}", paf_output);

    Ok(())
}
```

## Configuration

### Using the Builder Pattern

```rust
use fastga_rs::Config;

let config = Config::builder()
    .min_alignment_length(150)  // Minimum 150bp alignments
    .min_identity(0.8)          // 80% identity threshold
    .num_threads(8)             // Use 8 threads
    .soft_masking(true)         // Enable soft masking
    .build();
```

### Preset Configurations

```rust
// High sensitivity for distant homologs
let config = Config::high_sensitivity();

// Fast mode for closely related genomes
let config = Config::fast();

// Optimized for repetitive genomes
let config = Config::repetitive_genomes();
```

## Extended CIGAR Format

FastGA-RS produces alignments with extended CIGAR strings using:
- `=` for matches
- `X` for mismatches
- `I` for insertions
- `D` for deletions

This detailed format is essential for accurate variant calling and provides more information than traditional CIGAR strings that use `M` for both matches and mismatches.

## Advanced Usage

### Streaming API

The streaming API provides fine-grained control over alignment processing, allowing you to filter and process alignments as they're generated without storing them in memory or writing intermediate files.

#### Simple Streaming

```rust
use fastga_rs::streaming::align_streaming_simple;

let mut high_quality = Vec::new();

align_streaming_simple(genome1, genome2, |alignment| {
    if alignment.identity() > 0.95 {
        high_quality.push(alignment);
        true  // Keep processing
    } else {
        false  // Skip this alignment
    }
})?;
```

#### Advanced Filtering with StreamingAligner

```rust
use fastga_rs::streaming::StreamingAligner;
use fastga_rs::Config;

let mut aligner = StreamingAligner::new(Config::default());

// Compose multiple filters
aligner
    .filter_min_identity(0.9)        // Minimum 90% identity
    .filter_min_length(500)           // Minimum 500bp alignments
    .filter_query(|name| name.starts_with("chr"))  // Only chromosomes
    .filter_target(|name| !name.contains("_alt")); // Exclude alternates

let stats = aligner.align_files(genome1, genome2, |alignment| {
    // Process each alignment that passes all filters
    println!("{} -> {}: {:.1}%",
        alignment.query_name,
        alignment.target_name,
        alignment.identity() * 100.0);
    true
})?;

println!("Processed {} alignments, kept {}",
    stats.total_alignments,
    stats.kept_alignments);
```

#### Query-vs-All Mode

Efficiently align a single query against all targets:

```rust
let mut aligner = StreamingAligner::new(Config::default());

let hits = aligner.align_query_vs_all(
    query_sequence,
    target_database,
    |alignment| alignment.identity() > 0.85  // Keep high-identity hits
)?;

// Results are filtered and collected
for hit in &hits {
    println!("Hit: {} ({:.1}% identity)", hit.target_name, hit.identity() * 100.0);
}
```

#### Best Hit Per Query

Keep only the best alignment for each query:

```rust
use fastga_rs::streaming::BestHitFilter;

let mut best_filter = BestHitFilter::new();

align_streaming_simple(genome1, genome2, |alignment| {
    best_filter.process(alignment);
    true
})?;

let best_hits = best_filter.into_alignments();
```

#### Statistics Aggregation

Collect statistics without storing alignments:

```rust
let mut aligner = StreamingAligner::new(Config::default());

let mut total_bases = 0;
let mut identity_sum = 0.0;

aligner.aggregate(|alignment| {
    total_bases += alignment.query_end - alignment.query_start;
    identity_sum += alignment.identity();
});

let stats = aligner.align_files(genome1, genome2, |_| true)?;

let avg_identity = identity_sum / stats.total_alignments as f64;
println!("Average identity: {:.2}%", avg_identity * 100.0);
```

#### Stream to File

Stream alignments directly to a PAF file:

```rust
use fastga_rs::streaming::stream_to_paf;
use std::fs::File;
use std::io::BufWriter;

let file = File::create("alignments.paf")?;
let writer = BufWriter::new(file);

stream_to_paf(genome1, genome2, Config::default(), writer)?;
```

### Performance Benefits

The streaming API provides several performance advantages:

1. **Memory Efficiency**: Process TB-scale alignments without storing them
2. **Reduced I/O**: No intermediate .1aln files
3. **Early Termination**: Stop processing when you've found what you need
4. **Pipeline Integration**: Direct integration with downstream tools
5. **Parallel Processing**: Process alignments as they're generated

### Use Cases

#### Large-scale Genome Comparison
```rust
// Process human genome alignments without storing 100GB+ of alignments
stream_to_paf(human_genome1, human_genome2, config, output_file)?;
```

#### Contamination Screening
```rust
// Find and remove contaminant sequences
aligner.filter_target(|name| !known_contaminants.contains(name));
```

#### Rapid Database Search
```rust
// Stop after finding first high-quality hit
let mut found = false;
align_streaming_simple(query, database, |alignment| {
    if alignment.identity() > 0.99 {
        println!("Found perfect match: {}", alignment.target_name);
        found = true;
        false  // Stop processing
    } else {
        true  // Continue
    }
})?;
```

## Architecture

The library consists of several key components:

- **FFI Layer**: Low-level bindings to FastGA C code
- **Safe Wrappers**: Rust-safe abstractions over FFI
- **Configuration**: Flexible parameter management
- **Format Conversion**: PAF, PSL, and trace point formats
- **Streaming Interface**: Fine-grained alignment control

### Data Flow

```
FASTA files → GDB (genome database) → GIX (genome index) → Alignments → PAF/PSL output
```

## Performance Considerations

- **Thread Count**: Defaults to number of CPU cores
- **Memory Usage**: Scales with genome size and number of alignments
- **Temporary Files**: Uses system temp directory by default (configurable)
- **Compilation**: FastGA C code is compiled with `-O3` optimization

## Development

### Building from Source

```bash
git clone https://github.com/yourusername/fastga-rs
cd fastga-rs
cargo build --release
```

### Running Tests

```bash
cargo test
```

### Running Examples

```bash
cargo run --example basic_alignment
```

## Integration with Other Tools

FastGA-RS is designed to integrate well with other bioinformatics tools. The PAF output with extended CIGAR can be directly used with:
- Variant callers
- Alignment visualization tools
- Downstream analysis pipelines

## License

This project is dual-licensed under MIT OR Apache-2.0.

## Acknowledgments

FastGA was created by Gene Myers. This Rust wrapper aims to make FastGA's powerful alignment capabilities easily accessible to the Rust ecosystem.

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## Future Roadmap

- [ ] Complete FFI bindings for streaming alignment processing
- [ ] Implement query-vs-all optimization
- [ ] Add support for custom scoring matrices
- [ ] Integrate with Rust bioinformatics ecosystem (rust-bio, noodles)
- [ ] Add benchmarking suite
- [ ] Support for additional output formats