//! Example showing how to integrate FastGA with SweepGA to avoid FFI hanging issues.
//!
//! This uses the subprocess-based API which isolates each utility call in a separate process.

use anyhow::Result;
use fastga_rs::api::FastGA;
use fastga_rs::Config;
use std::path::Path;

fn main() -> Result<()> {
    // Configure FastGA
    let config = Config::builder()
        .num_threads(4)
        .min_alignment_length(100)
        .min_identity(0.7)
        .build();

    // Create the subprocess-based aligner (recommended for SweepGA)
    let aligner = FastGA::new(config)?;

    // Example paths (replace with actual files)
    let query_path = Path::new("query.fasta");
    let target_path = Path::new("target.fasta");

    // Check if files exist for demo
    if !query_path.exists() || !target_path.exists() {
        eprintln!("Demo files not found. This example shows the API usage.");
        eprintln!("To run with real data, provide query.fasta and target.fasta files.");
        return Ok(());
    }

    println!("Starting FastGA alignment (subprocess-based implementation)...");

    // Perform alignment - this will NOT hang unlike the FFI version
    match aligner.align_files(query_path, target_path) {
        Ok(alignments) => {
            println!("Alignment completed successfully!");
            println!("Found {} alignments", alignments.len());

            // Convert to PAF format for output
            if let Ok(paf) = alignments.to_paf() {
                println!("PAF output:");
                for line in paf.lines().take(5) {
                    println!("{}", line);
                }
                if alignments.len() > 5 {
                    println!("... ({} more alignments)", alignments.len() - 5);
                }
            }
        }
        Err(e) => {
            eprintln!("Alignment failed: {}", e);
            return Err(e.into());
        }
    }

    Ok(())
}

// Example of how to use this in SweepGA's sweep_algo.rs:
//
// ```rust
// use fastga_rs::api::FastGA;  // Use subprocess API instead of regular FastGA
// use fastga_rs::Config;
//
// pub fn run_fastga_alignment(
//     queries_path: &Path,
//     targets_path: &Path,
//     num_threads: usize,
// ) -> Result<Vec<AlignmentResult>> {
//     // Configure FastGA
//     let config = Config::builder()
//         .num_threads(num_threads as i32)
//         .min_alignment_length(100)
//         .min_identity(0.7)
//         .build();
//
//     // Use subprocess-based API to avoid hanging
//     let aligner = FastGA::new(config)?;
//
//     // This will NOT hang
//     let alignments = aligner.align_files(queries_path, targets_path)?;
//
//     // Process alignments...
//     Ok(process_alignments(alignments))
// }
// ```
