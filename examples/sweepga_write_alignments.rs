//! Example showing how to write alignments to disk for SweepGA integration.
//!
//! This demonstrates all the different ways to get alignment data to disk.

use anyhow::Result;
use fastga_rs::{Alignments, Config, FastGA};
use std::path::Path;

fn main() -> Result<()> {
    let query_path = Path::new("query.fasta");
    let target_path = Path::new("target.fasta");

    // Configure FastGA
    let config = Config::builder()
        .num_threads(4)
        .min_alignment_length(100)
        .min_identity(0.7)
        .build();

    let aligner = FastGA::new(config)?;

    println!("FastGA Alignment Output Options for SweepGA:\n");

    // Option 1: Align and write PAF directly
    println!("1. Direct PAF output:");
    let alignments = aligner.align_files(query_path, target_path)?;
    alignments.write_paf("output.paf")?;
    println!("   ✓ Wrote {} alignments to output.paf", alignments.len());

    // Option 2: Align to memory, then write various formats
    println!("\n2. In-memory processing with multiple output formats:");
    let alignments = aligner.align_files(query_path, target_path)?;

    // Write PAF format
    alignments.write_paf("output.paf")?;
    println!("   ✓ Wrote PAF format");

    // Write TSV for spreadsheet analysis
    alignments.write_tsv("output.tsv")?;
    println!("   ✓ Wrote TSV format (Excel-compatible)");

    // Write JSON for programmatic parsing
    alignments.write_json("output.json")?;
    println!("   ✓ Wrote JSON format");

    // Option 3: Process alignments before writing
    println!("\n3. Filtered output (only high-quality alignments):");
    let mut filtered = Alignments::new();
    filtered.alignments = alignments
        .alignments
        .iter()
        .filter(|a| a.identity() > 0.9)
        .cloned()
        .collect();
    filtered.write_paf("output_filtered.paf")?;
    println!("   ✓ Wrote {} high-identity alignments", filtered.len());

    // Option 4: Group by query and write separately
    println!("\n4. Per-query output files:");
    let groups = alignments.group_by_query();
    for (query_name, query_alignments) in groups.iter().take(3) {
        let filename = format!("{}_alignments.paf", query_name);
        let mut query_aligns = Alignments::new();
        query_aligns.alignments = query_alignments.iter().map(|&a| a.clone()).collect();
        query_aligns.write_paf(&filename)?;
        println!(
            "   ✓ Wrote {} alignments for query '{}'",
            query_alignments.len(),
            query_name
        );
    }

    // Option 5: Get summary statistics
    println!("\n5. Alignment summary:");
    let summary = alignments.summary();
    println!("   Total alignments: {}", summary.num_alignments);
    println!("   Mean identity: {:.1}%", summary.mean_identity * 100.0);
    println!(
        "   Identity range: {:.1}% - {:.1}%",
        summary.min_identity * 100.0,
        summary.max_identity * 100.0
    );

    println!("\n✅ All output methods demonstrated successfully!");

    Ok(())
}

// For per-query iteration (memory-efficient for large datasets):
#[allow(dead_code)]
fn process_queries_iteratively() -> Result<()> {
    use fastga_rs::QueryAlignmentIterator;

    let config = Config::default();
    let iter = QueryAlignmentIterator::new(
        Path::new("queries.fasta"),
        Path::new("targets.fasta"),
        config,
        1, // Buffer size
    )?;

    for query_result in iter {
        let query_set = query_result?;

        // Process each query's alignments individually
        println!(
            "Query '{}' has {} alignments",
            query_set.query_name,
            query_set.alignment_count()
        );

        // Write this query's alignments
        let mut aligns = Alignments::new();
        aligns.alignments = query_set.alignments;
        aligns.write_paf(format!("{}.paf", query_set.query_name))?;
    }

    Ok(())
}
