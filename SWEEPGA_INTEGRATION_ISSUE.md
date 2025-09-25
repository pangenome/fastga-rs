# FastGA-RS Integration Issue Report

## Problem Description

When integrating fastga-rs into sweepga, the FFI calls hang indefinitely. The library successfully compiles and embeds FastGA (confirmed by `fastga_main` symbol in binary), but execution never completes.

## How to Reproduce

1. **In sweepga directory**, the integration is set up as:

```rust
// Cargo.toml
[dependencies]
fastga-rs = { git = "https://github.com/pangenome/fastga-rs.git", rev = "82b6b186bb6a8745a14c352f2350eaa2e3a8650c" }
```

```rust
// src/fastga_integration.rs
use fastga_rs::{FastGA, Config};

let config = Config::builder()
    .num_threads(4)
    .build();

let aligner = FastGA::new(config)?;
let alignments = aligner.align_files(queries_path, targets_path)?;  // <- HANGS HERE
```

2. **Test with minimal data**:
```bash
# Create test file
echo -e ">seq1\nACGTACGTACGT" > test.fa

# Try to run (will hang)
cd ~/sweepga
cargo run -- test.fa -t 1

# Process hangs at "Running FastGA alignment..." and must be killed
```

3. **Even simpler test directly in fastga-rs**:
```rust
// Add this test to fastga-rs/tests/integration_test.rs
#[test]
fn test_minimal_alignment() {
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    let dir = tempdir().unwrap();
    let test_fa = dir.path().join("test.fa");

    let mut file = File::create(&test_fa).unwrap();
    writeln!(file, ">seq1\nACGTACGT").unwrap();

    let config = Config::builder()
        .num_threads(1)
        .build();

    let aligner = FastGA::new(config).unwrap();

    println!("Starting alignment...");
    let result = aligner.align_files(&test_fa, &test_fa);  // Self-alignment
    println!("Alignment complete!");

    assert!(result.is_ok());
}
```

Run with: `cargo test test_minimal_alignment -- --nocapture`

## Observed Behavior

- The `FastGA::new()` call succeeds
- The `align_files()` call never returns
- No error messages are produced
- Process CPU usage goes to 0% (not spinning, appears to be waiting/blocked)
- Must be terminated with Ctrl-C or timeout

## Suspected Issues

1. **Temp directory/file handling**: The orchestrator may be having issues creating or accessing temp files/directories

2. **FFI blocking**: The C functions called via FFI might be:
   - Waiting for input that never comes
   - Deadlocked on a mutex/resource
   - Having path resolution issues

3. **Missing initialization**: Some FastGA initialization step might be missing when called via FFI vs command line

## Suggested API Improvements

### 1. Add timeout/async support
```rust
// Allow timeout
let alignments = aligner.align_files_timeout(
    queries,
    targets,
    Duration::from_secs(30)
)?;

// Or async version
let alignments = aligner.align_files_async(queries, targets).await?;
```

### 2. Add progress callback for debugging
```rust
let aligner = FastGA::new(config)?
    .on_progress(|stage, message| {
        eprintln!("FastGA: {} - {}", stage, message);
    });
```

### 3. Expose intermediate steps for debugging
```rust
// Instead of one monolithic align_files(), expose the pipeline:
let aligner = FastGA::new(config)?;
let query_db = aligner.prepare_database(queries)?;  // Test this works
let target_db = aligner.prepare_database(targets)?; // Test this works
let alignments = aligner.align_databases(&query_db, &target_db)?; // Test this
```

### 4. Add dry-run or validate mode
```rust
// Check if alignment would work without actually running
let validation = aligner.validate_alignment(queries, targets)?;
if validation.is_ok() {
    let alignments = aligner.align_files(queries, targets)?;
}
```

## Quick Debug Steps

1. Add logging to `orchestrator.rs` at each step:
```rust
eprintln!("Creating temp dir...");
eprintln!("Preparing query database...");
eprintln!("Preparing target database...");
eprintln!("Running FastGA main...");
```

2. Check if the issue is with file paths - try absolute vs relative paths

3. Verify temp directory permissions and that files are actually being created

4. Try running the embedded `fastga_main` directly with simple test to confirm FFI works:
```rust
// Minimal FFI test
let args = vec!["fastga", "--version"];
let result = unsafe { fastga_main(args.len() as i32, args.as_ptr()) };
println!("Result: {}", result);
```

## Contact

If you need help testing or want to see the exact integration code, check:
- `/home/erik/sweepga/src/fastga_integration.rs`
- `/home/erik/sweepga/src/main.rs` (lines 187-236 for FastGA mode handling)

The integration is clean and minimal - just Config with threads, then `align_files()`. The hanging happens immediately on the align_files call.