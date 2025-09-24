//! Example demonstrating basic FastGA alignment usage.

use fastga_rs::{FastGA, Config};
use std::path::Path;
use anyhow::Result;

fn main() -> Result<()> {
    // Create aligner with default configuration
    println!("Creating FastGA aligner with default configuration...");
    let aligner = FastGA::new(Config::default())?;

    // Example 1: Align two genome files
    println!("\n=== Example 1: Basic file alignment ===");
    let genome1 = Path::new("data/genome1.fasta");
    let genome2 = Path::new("data/genome2.fasta");

    // Note: In a real scenario, these files would exist
    // let alignments = aligner.align_files(genome1, genome2)?;
    // println!("Found {} alignments", alignments.len());

    // Example 2: Using custom configuration
    println!("\n=== Example 2: Custom configuration ===");
    let config = Config::builder()
        .min_alignment_length(200)  // Minimum 200bp alignments
        .min_identity(0.85)         // 85% identity threshold
        .num_threads(4)             // Use 4 threads
        .build();

    let custom_aligner = FastGA::new(config)?;
    println!("Created aligner with custom settings");

    // Example 3: Using preset configurations
    println!("\n=== Example 3: Preset configurations ===");

    // High sensitivity for distant homologs
    let sensitive_aligner = FastGA::new(Config::high_sensitivity())?;
    println!("Created high-sensitivity aligner for distant homologs");

    // Fast mode for closely related genomes
    let fast_aligner = FastGA::new(Config::fast())?;
    println!("Created fast aligner for closely related genomes");

    // For repetitive genomes
    let repeat_aligner = FastGA::new(Config::repetitive_genomes())?;
    println!("Created aligner optimized for repetitive genomes");

    // Example 4: Align sequences from memory
    println!("\n=== Example 4: Aligning sequences from memory ===");
    let seq1 = b">seq1\nACGTACGTACGTACGT\n";
    let seq2 = b">seq2\nACGTACGTACGTACGT\n";

    // let alignments = aligner.align_sequences(seq1, seq2)?;
    // println!("Aligned sequences from memory");

    // Example 5: Processing alignments
    println!("\n=== Example 5: Processing alignment results ===");

    // Example of how to process alignments when available:
    /*
    for alignment in alignments.iter() {
        println!("Query: {} -> Target: {}",
            alignment.query_name,
            alignment.target_name);
        println!("  Identity: {:.1}%", alignment.identity() * 100.0);
        println!("  Query coverage: {:.1}%", alignment.query_coverage() * 100.0);
        println!("  CIGAR: {}", alignment.cigar);
    }

    // Filter high-quality alignments
    let mut filtered = alignments;
    filtered.filter(|a| a.identity() > 0.9 && a.query_coverage() > 0.5);
    println!("High-quality alignments: {}", filtered.len());

    // Convert to PAF format
    let paf_output = filtered.to_paf()?;
    println!("PAF output:\n{}", paf_output);
    */

    // Example 6: Streaming API (when implemented)
    println!("\n=== Example 6: Streaming API (future) ===");
    println!("The streaming API will allow processing alignments as they're generated:");
    println!("- Filter alignments on-the-fly");
    println!("- Reduce memory usage for large-scale alignments");
    println!("- Implement custom processing pipelines");

    /*
    // Future streaming API example:
    let mut count = 0;
    aligner.align_streaming(genome1, genome2, |alignment| {
        if alignment.identity() > 0.95 {
            count += 1;
            // Process high-identity alignment
        }
        true  // Continue processing
    })?;
    println!("Found {} high-identity alignments", count);
    */

    Ok(())
}