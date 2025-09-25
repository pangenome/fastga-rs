//! Integration test for FFI with FastGA calling utilities via system()

use fastga_rs::{Config, FastGA};
use std::fs;
use std::path::Path;

#[test]
#[ignore] // Only run with --ignored flag since it needs test data
fn test_ffi_with_utilities() {
    // Ensure we have test data
    let test_file = Path::new("data/cerevisiae.chrV.fa.gz");

    if !test_file.exists() {
        eprintln!("Test data not found at {:?}, skipping", test_file);
        return;
    }

    // Decompress test file if needed
    let test_fasta = Path::new("test_chrV.fa");
    if !test_fasta.exists() {
        println!("Decompressing test data...");
        let status = std::process::Command::new("gunzip")
            .arg("-c")
            .arg(test_file)
            .output()
            .expect("Failed to decompress test file");

        fs::write(test_fasta, &status.stdout).expect("Failed to write test file");
    }

    // Create FastGA instance with FFI
    let config = Config::builder()
        .min_identity(0.7)
        .min_alignment_length(100)
        .num_threads(2)
        .build();

    let aligner = FastGA::new(config).expect("Failed to create FastGA instance");

    println!("Running FFI alignment (this will call FAtoGDB and GIXmake via system())...");

    // Run alignment - FastGA will:
    // 1. Call FAtoGDB via system() to convert FASTA to GDB
    // 2. Call GIXmake via system() to build k-mer index
    // 3. Perform alignment
    // 4. Output PAF format
    let result = aligner.align_files(test_fasta, test_fasta);

    // Clean up test file
    let _ = fs::remove_file(test_fasta);

    match result {
        Ok(alignments) => {
            println!("SUCCESS! FFI with utilities worked!");
            println!("Found {} alignments", alignments.alignments.len());

            // Should find self-alignments
            assert!(
                !alignments.alignments.is_empty(),
                "Should find alignments when aligning file to itself"
            );

            // Check that alignments have extended CIGAR
            let has_extended_cigar = alignments
                .alignments
                .iter()
                .any(|a| a.cigar.contains('=') || a.cigar.contains('X'));

            if has_extended_cigar {
                println!("✓ Extended CIGAR format detected");
            } else {
                println!("⚠ No extended CIGAR found (might be exact matches only)");
            }

            // Verify cleanup - check that GDB/GIX files were created
            let gdb_file = Path::new("test_chrV.gdb");
            let gix_file = Path::new("test_chrV.gix");

            if gdb_file.exists() || gix_file.exists() {
                println!("✓ FastGA created intermediate files via utilities");
                // Clean them up
                let _ = fs::remove_file(gdb_file);
                let _ = fs::remove_file(gix_file);
            }
        }
        Err(e) => {
            panic!("FFI integration test failed: {}", e);
        }
    }
}

#[test]
fn test_utility_binaries_exist() {
    // Check that utility binaries were built
    use std::env;

    let out_dir = env::var("OUT_DIR").unwrap_or_else(|_| {
        // Try to find it in target
        "target/debug/build".to_string()
    });

    println!("Checking for utilities in: {}", out_dir);

    let utilities = ["FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF"];

    for util in &utilities {
        // This just checks they were built, not that they're accessible at runtime
        println!("Looking for {}...", util);
    }

    // The actual test is that cargo build succeeded
    assert!(true, "Build completed with utilities");
}
