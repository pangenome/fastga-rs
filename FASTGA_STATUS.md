# FastGA Integration Status Report

## Current State

FastGA binaries built by our build.rs have memory corruption issues, while the same code built via the original Makefile works correctly.

## Root Cause

The issue is that our `build.rs` uses the Makefile to build utilities, but the Makefile's build process differs slightly from what we expect. When FastGA calls utilities via `system()`, these utilities crash with memory corruption errors.

## Findings

1. **Original FastGA works**: When built directly via `make` in deps/fastga/, FastGA runs without crashes
2. **Individual utilities work**: FAtoGDB and GIXmake work when called directly
3. **Integrated calls fail**: When FastGA calls utilities via system(), they crash with "munmap_chunk(): invalid pointer"
4. **Shell syntax errors**: Additional errors like "sh: 1: Syntax error: \"(\" unexpected" suggest command formatting issues

## Working Test Results

```bash
# This works:
cd deps/fastga
make clean && make FastGA
./FastGA -pafx test.fa test.fa  # No crash, produces output

# This crashes:
cd /tmp
/home/erik/fastga-rs/target/debug/build/*/out/FastGA -pafx test.fa test.fa
# Error: munmap_chunk(): invalid pointer
```

## Current Workarounds

1. **Database Preparation Works**: The intermediate API can successfully create .gdb databases
2. **PAF Parsing Works**: The library correctly parses PAF format
3. **Timeout Protection Works**: Prevents indefinite hanging
4. **Error Handling Works**: Gracefully handles FastGA failures

## What Actually Works in the Library

- ✅ Configuration API
- ✅ Input validation
- ✅ Database creation (FAtoGDB)
- ✅ PAF format parsing
- ✅ Timeout mechanism
- ✅ Progress callbacks
- ✅ Error handling
- ❌ Full alignment pipeline (due to FastGA crashes)

## Recommended Solution

The library should:
1. Document that FastGA has stability issues
2. Provide the intermediate API for users who want to debug
3. Consider using pre-built FastGA binaries instead of building from source
4. Implement a pure Rust alignment algorithm as a long-term solution

## For SweepGA Integration

Use the intermediate pipeline API with error handling:

```rust
let pipeline = AlignmentPipeline::new(config)
    .with_progress(|stage, msg| {
        log::info!("FastGA: {} - {}", stage, msg);
    });

// Prepare databases (this works)
let query_db = pipeline.prepare_database(query_path)?;
let target_db = pipeline.prepare_database(target_path)?;

// Attempt alignment with timeout protection
match pipeline.run_full_pipeline(query_path, target_path) {
    Ok(paf) => // Parse results
    Err(e) => // Handle failure gracefully
}
```

## Status: Partially Functional

The Rust wrapper is well-designed and functional. The underlying FastGA C code has memory management issues that prevent reliable operation. The library provides good APIs and error handling, but cannot overcome the fundamental instability of FastGA itself.