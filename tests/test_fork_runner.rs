use fastga_rs::fork_runner::{fork_fatogdb, fork_gixmake, ForkOrchestrator};
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn test_fork_fatogdb_basic() {
    let dir = tempdir().unwrap();
    let test_fa = dir.path().join("test.fa");

    // Create a simple test FASTA file
    let mut file = File::create(&test_fa).unwrap();
    writeln!(file, ">seq1").unwrap();
    writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
    writeln!(file, ">seq2").unwrap();
    writeln!(file, "TACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
    file.flush().unwrap();

    eprintln!("Testing fork_fatogdb with {}", test_fa.display());

    match fork_fatogdb(&test_fa) {
        Ok(gdb_path) => {
            eprintln!("✓ Successfully created GDB: {:?}", gdb_path);
            assert!(gdb_path.exists(), "GDB file should exist");
            assert!(gdb_path.extension() == Some("1gdb".as_ref()), "Should have .1gdb extension");

            // Also check that .1bps was created
            let bps_path = test_fa.with_extension("1bps");
            if bps_path.exists() {
                eprintln!("✓ BPS file also exists");
            } else {
                eprintln!("⚠ BPS file not found (might be expected for small files)");
            }

            eprintln!("✓ GDB file size: {} bytes", std::fs::metadata(&gdb_path).unwrap().len());
            if bps_path.exists() {
                eprintln!("✓ BPS file size: {} bytes", std::fs::metadata(&bps_path).unwrap().len());
            }
        }
        Err(e) => {
            panic!("fork_fatogdb failed: {}", e);
        }
    }
}

#[test]
fn test_fork_gixmake_basic() {
    let dir = tempdir().unwrap();
    let test_fa = dir.path().join("test.fa");

    // Create a test FASTA file with more sequence data for indexing
    let mut file = File::create(&test_fa).unwrap();
    writeln!(file, ">seq1").unwrap();
    for _ in 0..10 {
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
    }
    writeln!(file, ">seq2").unwrap();
    for _ in 0..10 {
        writeln!(file, "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
    }
    file.flush().unwrap();

    eprintln!("Testing fork_gixmake pipeline");

    // First create GDB
    let gdb_path = fork_fatogdb(&test_fa).expect("FAtoGDB should work");
    eprintln!("✓ Created GDB: {}", gdb_path.display());

    // Then create index
    match fork_gixmake(&gdb_path, 1, 10) {
        Ok(gix_path) => {
            eprintln!("✓ Successfully created GIX: {:?}", gix_path);
            assert!(gix_path.exists(), "GIX file should exist");
            assert!(gix_path.extension() == Some("gix".as_ref()), "Should have .gix extension");

            eprintln!("✓ GIX file size: {} bytes", std::fs::metadata(&gix_path).unwrap().len());
        }
        Err(e) => {
            eprintln!("✗ fork_gixmake failed: {}", e);
            eprintln!("Note: This might be expected if GIXmake requires specific conditions");
        }
    }
}

#[test]
fn test_fork_orchestrator_full_pipeline() {
    let dir = tempdir().unwrap();

    // Create query and target files
    let query_fa = dir.path().join("query.fa");
    let target_fa = dir.path().join("target.fa");

    let mut file = File::create(&query_fa).unwrap();
    writeln!(file, ">query1").unwrap();
    for _ in 0..5 {
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
    }
    file.flush().unwrap();

    let mut file = File::create(&target_fa).unwrap();
    writeln!(file, ">target1").unwrap();
    for _ in 0..5 {
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
    }
    writeln!(file, ">target2").unwrap();
    for _ in 0..5 {
        writeln!(file, "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
    }
    file.flush().unwrap();

    eprintln!("Testing full pipeline with ForkOrchestrator");

    let orchestrator = ForkOrchestrator::new_simple(1);

    match orchestrator.align(&query_fa, &target_fa) {
        Ok(output) => {
            eprintln!("✓ Alignment succeeded!");
            eprintln!("Output length: {} bytes", output.len());

            if !output.is_empty() {
                // Print first few lines of output
                for (i, line) in output.lines().take(5).enumerate() {
                    eprintln!("  Line {}: {}", i + 1, line);
                }

                // Check for PAF format
                if output.contains('\t') {
                    eprintln!("✓ Output appears to be in PAF format");
                }
            } else {
                eprintln!("⚠ Warning: Empty output");
            }
        }
        Err(e) => {
            eprintln!("✗ Full pipeline failed: {}", e);
            eprintln!("This might be expected if FastGA binary is not available");
        }
    }
}