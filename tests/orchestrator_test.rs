//! Test the pure Rust orchestration of FastGA

use fastga_rs::{Config, FastGA};
use std::fs;
use std::path::Path;

#[test]
fn test_orchestrator_builds() {
    // Just test that we can create an instance
    let aligner = FastGA::new(Config::default());
    assert!(aligner.is_ok());
}

#[test]
#[ignore] // Run with --ignored to test with real data
fn test_orchestrator_alignment() {
    let test_file = Path::new("data/cerevisiae.chrV.fa.gz");

    if !test_file.exists() {
        eprintln!("Test data not found, skipping");
        return;
    }

    // Decompress for testing
    let test_fa = Path::new("test_orchestrator.fa");
    let output = std::process::Command::new("gunzip")
        .arg("-c")
        .arg(test_file)
        .output()
        .expect("Failed to decompress");

    fs::write(test_fa, &output.stdout).expect("Failed to write test file");

    // Create aligner with our orchestrator
    let config = Config::builder()
        .min_identity(0.7)
        .min_alignment_length(100)
        .num_threads(2)
        .build();

    let aligner = FastGA::new(config).expect("Failed to create aligner");

    println!("Testing orchestrator-based alignment...");
    println!("This calls fatogdb_main() and gixmake_main() directly via FFI");

    let result = aligner.align_files(test_fa, test_fa);

    // Clean up
    let _ = fs::remove_file(test_fa);
    let _ = fs::remove_file("test_orchestrator.gdb");
    let _ = fs::remove_file("test_orchestrator.gix");

    match result {
        Ok(alignments) => {
            println!("✓ Orchestrator SUCCESS!");
            println!("  Found {} alignments", alignments.alignments.len());
            println!("  All FFI calls worked correctly");

            assert!(!alignments.alignments.is_empty());
        }
        Err(e) => {
            panic!("Orchestrator failed: {}", e);
        }
    }
}
