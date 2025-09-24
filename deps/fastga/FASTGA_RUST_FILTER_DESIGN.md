# FASTGA + Rust Filtering Library Design Document

## Overview
Design for a standalone Rust library that combines FASTGA's fast alignment engine with wfmash's sophisticated filtering algorithms. This would be a completely separate system from wfmash, reimplementing the filtering logic in Rust.

## Architecture

### Core Components

1. **FASTGA Integration Layer**
   - Call FASTGA as external process or via FFI
   - Parse PAF output stream
   - Convert to internal compact mapping format

2. **Compact Mapping Structure (32 bytes)**
   ```rust
   #[repr(C, packed)]
   pub struct CompactMapping {
       ref_seq_id: u32,
       ref_start_pos: u32,
       query_start_pos: u32,
       block_length: u32,
       n_merged: u32,
       conserved_sketches: u32,
       nuc_identity: u16,      // 0-10000 for 0.00-100.00%
       flags: u8,               // bits: 0=strand, 1=discard, 2=overlapped
       kmer_complexity: u8,     // 0-100
   }
   ```

3. **Filtering Pipeline** (reimplemented from wfmash)
   - `merge_mappings_in_range()` - Chain mappings within gap distance
   - `filter_weak_mappings()` - Remove short merged blocks
   - `filter_by_group()` - Keep top N mappings per reference group
   - `filter_false_high_identity()` - Remove length-mismatched alignments
   - `sparsify_mappings()` - Reduce redundant mappings
   - `filter_by_scaffolds()` - Complex scaffold-aware filtering
   - `one_to_one_filter()` - Reciprocal best hit filtering

## Implementation Strategy

### Phase 1: FASTGA Wrapper
```rust
// fastga-filter/src/fastga.rs
pub struct FastGARunner {
    target_path: PathBuf,
    query_path: PathBuf,
    params: FastGAParams,
}

impl FastGARunner {
    pub fn run(&self) -> Result<PafReader> {
        // Execute FASTGA process
        // Return streaming PAF reader
    }
}
```

### Phase 2: Filtering Engine
```rust
// fastga-filter/src/filter.rs
pub struct FilterEngine {
    mappings: Vec<CompactMapping>,
    chain_info: Option<Vec<ChainInfo>>,
    aux_data: Option<Vec<MappingAuxData>>,
}

impl FilterEngine {
    pub fn from_paf(reader: PafReader) -> Self {
        // Parse PAF, convert to compact format
    }

    pub fn apply_filters(&mut self, config: FilterConfig) {
        // Apply filtering pipeline in order
    }

    pub fn to_paf(&self) -> String {
        // Convert back to PAF format
    }
}
```

### Phase 3: Comparison Framework
```rust
// fastga-filter/benches/compare.rs
fn benchmark_against_wfmash(test_data: &TestData) {
    // Run both systems
    // Compare outputs
    // Measure performance
}
```

## Key Design Decisions

### Memory Efficiency
- **32-byte compact struct** - Matches wfmash's memory layout exactly
- **Separate auxiliary data** - Chain info only allocated when needed
- **Bit-packed flags** - Multiple booleans in single byte
- **Scaled values** - Identity as uint16 (0-10000), complexity as uint8 (0-100)

### Filtering Algorithms (from wfmash)
1. **Chain Merging**: Merge mappings within specified gap distance
2. **Group Filtering**: Keep top N mappings per reference sequence/group
3. **Scaffold Filtering**: Complex filtering considering scaffold structure
4. **One-to-One**: Reciprocal best hit filtering for 1:1 ortholog detection

### Integration Options

#### Option A: Process-based
- Run FASTGA as subprocess
- Stream PAF output through pipe
- Parse and filter in Rust
- Output filtered PAF

#### Option B: FFI-based
- Link FASTGA as C library
- Direct memory access to alignment structures
- Convert Overlap structs to CompactMapping
- More efficient but more complex

## Advantages Over Current Approaches

1. **Clean Separation**: Completely independent from wfmash codebase
2. **Memory Safety**: Rust's ownership model prevents memory errors
3. **Modern Tooling**: Cargo for deps, built-in testing/benchmarking
4. **Performance**: Zero-cost abstractions, easy parallelization with Rayon
5. **Maintainability**: Single language, clear interfaces

## Comparison Strategy

### Test Suite
1. Use same test genomes (aEleCoq1 vs rDibSmi1)
2. Run three configurations:
   - FASTGA alone (baseline)
   - FASTGA + Rust filtering
   - wfmash (for comparison)

### Metrics
- Memory usage (peak RSS)
- Runtime (wall clock)
- Output similarity (mapping counts, positions)
- Filtering effectiveness (reduction ratios)

## File Structure
```
fastga-filter/
├── Cargo.toml
├── src/
│   ├── lib.rs           # Public API
│   ├── fastga.rs        # FASTGA runner
│   ├── parser.rs        # PAF parser
│   ├── filter.rs        # Filtering algorithms
│   ├── compact.rs       # Compact mapping types
│   └── output.rs        # Output formatting
├── benches/
│   └── compare.rs       # Comparison benchmarks
└── tests/
    └── integration.rs   # Integration tests
```

## Build Integration
```toml
# Cargo.toml
[dependencies]
nom = "7.1"              # PAF parsing
rayon = "1.7"            # Parallelization
kiddo = "2.0"            # KD-tree for spatial ops
memmap2 = "0.5"          # Memory-mapped file I/O
anyhow = "1.0"           # Error handling

[build-dependencies]
cc = "1.0"               # For building FASTGA if using FFI

[dev-dependencies]
criterion = "0.5"        # Benchmarking
proptest = "1.0"         # Property testing
```

## Usage Example
```rust
use fastga_filter::{FastGARunner, FilterEngine, FilterConfig};

fn main() -> Result<()> {
    // Run FASTGA
    let runner = FastGARunner::new("target.fa", "query.fa")
        .min_identity(0.85)
        .min_length(1000);

    let paf_reader = runner.run()?;

    // Apply filters
    let mut engine = FilterEngine::from_paf(paf_reader)?;

    let config = FilterConfig::default()
        .chain_gap(10000)
        .max_mappings_per_group(100)
        .enable_scaffold_filter(5000);

    engine.apply_filters(config);

    // Output filtered results
    println!("{}", engine.to_paf());

    Ok(())
}
```

## Next Steps
1. Create repository structure
2. Implement PAF parser
3. Port filtering algorithms from wfmash
4. Create benchmark suite
5. Optimize and compare performance