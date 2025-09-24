/// Example demonstrating the self-contained FastGA integration.
/// This example shows that fastga-rs embeds the FastGA binary and requires
/// no external dependencies - just a single Rust binary!

use anyhow::Result;
use fastga_rs::{FastGA, Config};
use std::path::Path;

fn main() -> Result<()> {
    println!("FastGA-RS Self-Contained Binary Example");
    println!("========================================");
    println!();
    println!("This example demonstrates that fastga-rs is completely self-contained.");
    println!("No external FastGA binary is required - everything is embedded!");
    println!();

    // Create temporary test sequences
    let seq1 = b">Sequence1\nACGTACGTACGTACGTACGT\n";
    let seq2 = b">Sequence2\nACGTACGTACGTACGTACGT\n";

    // Create aligner with default config
    let aligner = FastGA::new(Config::default())?;
    println!("✓ Created FastGA aligner (using embedded binary)");

    // Align sequences
    println!("✓ Running alignment...");
    let alignments = aligner.align_sequences(seq1, seq2)?;

    // Show results
    println!("✓ Alignment complete!");
    println!();
    println!("Alignments found: {}", alignments.alignments.len());

    if let Some(first) = alignments.alignments.first() {
        println!("First alignment:");
        println!("  Query: {}", first.query_name);
        println!("  Target: {}", first.target_name);
        println!("  Identity: {:.2}%", first.identity() * 100.0);
        println!("  CIGAR: {}", &first.cigar[..50.min(first.cigar.len())]);
        if first.cigar.len() > 50 {
            println!("         ...");
        }
    }

    println!();
    println!("SUCCESS: FastGA-RS is working with embedded binaries!");
    println!("This is a single self-contained Rust binary with no external dependencies.");

    Ok(())
}