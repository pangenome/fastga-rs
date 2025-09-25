# FastGA-RS Integration Fix Documentation

## Issue Resolution

The hanging issue reported in `SWEEPGA_INTEGRATION_ISSUE.md` has been addressed. The root cause was that the FFI calls to FastGA's main functions were blocking indefinitely. Additionally, there were memory corruption issues when running the compiled FastGA binaries.

## Solution Implemented

We've implemented a multi-layered solution with fallback mechanisms:

### 1. **Intermediate Pipeline API** (Recommended)

The new `intermediate::AlignmentPipeline` provides step-by-step control over the alignment process:

```rust
use fastga_rs::intermediate::AlignmentPipeline;
use fastga_rs::Config;

let config = Config::builder()
    .num_threads(4)
    .min_alignment_length(100)
    .build();

let pipeline = AlignmentPipeline::new(config)
    .with_progress(|stage, msg| {
        eprintln!("[Progress] {}: {}", stage, msg);
    });

// Option 1: Run full pipeline
let paf_output = pipeline.run_full_pipeline(query_path, target_path)?;

// Option 2: Run individual steps for debugging
pipeline.validate_inputs(query_path, target_path)?;
let query_db = pipeline.prepare_database(query_path)?;
let target_db = pipeline.prepare_database(target_path)?;
pipeline.create_index(&query_db)?;
pipeline.create_index(&target_db)?;
let paf_output = pipeline.align_databases(&query_db, &target_db)?;
```

### 2. **Timeout Support**

Prevent indefinite hanging with the new timeout API:

```rust
use fastga_rs::timeout::TimeoutAligner;
use std::time::Duration;

let aligner = TimeoutAligner::new(config)
    .with_timeout(Duration::from_secs(30))
    .with_progress(|stage, msg| {
        eprintln!("[FastGA] {}: {}", stage, msg);
    });

let alignments = aligner.align_files(queries, targets)?;
```

### 3. **Simple Runner Fallback**

If FFI continues to have issues, use the process-based runner:

```rust
use fastga_rs::simple_runner::run_fastga_simple;

let paf_output = run_fastga_simple(
    query_path,
    target_path,
    num_threads,
    min_length,
    min_identity,
)?;
```

## Known Issues and Workarounds

### Issue 1: Memory Corruption
**Symptom**: "munmap_chunk(): invalid pointer" error
**Cause**: FastGA's C code has memory management issues with certain inputs
**Workaround**: Use well-formed FASTA files with sequences longer than 30bp

### Issue 2: FFI Hanging
**Symptom**: FFI calls never return
**Cause**: FastGA's main functions expect specific initialization
**Workaround**: Use the intermediate pipeline or simple runner instead of direct FFI

### Issue 3: Missing Utilities
**Symptom**: "sh: 1: Syntax error" when FastGA calls utilities
**Cause**: FastGA uses `system()` to call FAtoGDB, GIXmake
**Workaround**: Ensure all utilities are in the same directory as FastGA

## Integration Example for SweepGA

```rust
// In sweepga/src/fastga_integration.rs

use fastga_rs::intermediate::AlignmentPipeline;
use fastga_rs::Config;
use std::path::Path;
use std::time::Duration;

pub fn run_alignment(
    queries: &Path,
    targets: &Path,
    num_threads: usize,
) -> anyhow::Result<Vec<Alignment>> {
    // Create config
    let config = Config::builder()
        .num_threads(num_threads)
        .min_alignment_length(100)
        .min_identity(0.85)
        .build();

    // Use pipeline with progress reporting
    let pipeline = AlignmentPipeline::new(config)
        .with_progress(|stage, msg| {
            log::info!("FastGA {}: {}", stage, msg);
        });

    // Run with error handling
    match pipeline.run_full_pipeline(queries, targets) {
        Ok(paf_output) => {
            // Parse PAF output into alignments
            fastga_rs::Alignments::from_paf(&paf_output)
                .map(|a| a.alignments)
        }
        Err(e) => {
            log::error!("FastGA alignment failed: {}", e);

            // Fallback to simple runner
            log::info!("Trying fallback method...");
            let paf = fastga_rs::simple_runner::run_fastga_simple(
                queries, targets, num_threads, 100, Some(0.85)
            )?;

            fastga_rs::Alignments::from_paf(&paf)
                .map(|a| a.alignments)
        }
    }
}
```

## Testing

Run the new test suite to verify the fixes work:

```bash
# Test intermediate pipeline
cargo test test_intermediate_pipeline -- --nocapture

# Test timeout functionality
cargo test test_timeout_aligner -- --nocapture

# Test with real data (if available)
cargo test test_with_real_data -- --ignored --nocapture
```

## Building from Source

The build process now correctly compiles FastGA and its utilities:

```bash
# Clean build
cargo clean
cargo build --release

# Verify utilities are built
ls -la target/release/build/*/out/
# Should see: FastGA, FAtoGDB, GIXmake, GIXrm, ALNtoPAF
```

## Environment Variables

- `CARGO_MANIFEST_DIR`: Set this to the fastga-rs directory if running from another location
- `OUT_DIR`: Build script sets this to the output directory containing utilities

## Next Steps

1. The FFI approach needs more work to properly initialize FastGA's internal state
2. Consider wrapping FastGA in a more robust C interface designed for FFI
3. Investigate the memory corruption issues in FastGA's C code
4. Add more comprehensive test coverage with various input types

## Contact

For issues, please report to https://github.com/pangenome/fastga-rs/issues with:
- The exact error message
- Input file characteristics (size, number of sequences)
- Output of `cargo test test_intermediate_pipeline -- --nocapture`