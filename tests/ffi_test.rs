//! Test that the FFI bindings work correctly

use fastga_rs::{Config, FastGA};
use std::path::Path;

#[test]
#[ignore] // Requires test data files
fn test_ffi_alignment() {
    // Create aligner with FFI backend
    let config = Config::builder()
        .min_identity(0.7)
        .min_alignment_length(100)
        .num_threads(2)
        .build();

    let aligner = FastGA::new(config).unwrap();

    // Test files
    let query = Path::new("data/cerevisiae.chrV.fa.gz");
    let target = Path::new("data/cerevisiae.chrV.fa.gz");

    if !query.exists() || !target.exists() {
        eprintln!("Test data not found, skipping test");
        return;
    }

    // Run alignment using FFI
    let result = aligner.align_files(query, target);

    match result {
        Ok(alignments) => {
            println!(
                "FFI Success! Found {} alignments",
                alignments.alignments.len()
            );
            assert!(
                !alignments.alignments.is_empty(),
                "Should find some alignments"
            );

            // Check first alignment has extended CIGAR
            if let Some(first) = alignments.alignments.first() {
                println!(
                    "First alignment: {} -> {}",
                    first.query_name, first.target_name
                );
                println!("CIGAR: {}", first.cigar);

                // Extended CIGAR should have '=' or 'X' operators
                assert!(
                    first.cigar.contains('=') || first.cigar.contains('X'),
                    "Should have extended CIGAR with = or X operators"
                );
            }
        }
        Err(e) => {
            eprintln!("FFI alignment failed: {}", e);
            panic!("FFI test failed");
        }
    }
}

#[test]
fn test_ffi_loads() {
    // Just test that we can create an aligner (library loads)
    let aligner = FastGA::new(Config::default());
    assert!(
        aligner.is_ok(),
        "Should be able to create FastGA instance with FFI"
    );
}
