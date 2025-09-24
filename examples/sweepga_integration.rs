//! Example of how to integrate FastGA-rs into SweepGA
//!
//! This shows how SweepGA could use FastGA-rs as a library to:
//! 1. Generate alignments on-demand
//! 2. Apply immediate plane sweep filtering
//! 3. Avoid writing massive intermediate PAF files

use fastga_rs::{Config, integrated::IntegratedAligner};
use fastga_rs::plane_sweep::PlaneSweepConfig;
use std::path::Path;

/// Example of how SweepGA's main processing loop could work
fn sweepga_main_with_integrated_aligner(
    query_file: &Path,
    target_file: &Path,
    output_file: &Path,
) -> Result<(), Box<dyn std::error::Error>> {

    // Configure FastGA for the alignment
    let align_config = Config::builder()
        .min_identity(0.85)
        .min_alignment_length(5000)
        .num_threads(8)
        .build();

    // Configure aggressive plane sweep for repetitive genomes
    let sweep_config = PlaneSweepConfig {
        max_per_query: 100,     // Limit alignments per query
        max_per_target: 1,      // 1:1 mapping
        min_identity: 0.90,
        min_length: 5000,
        max_overlap: 0.5,
    };

    // Create integrated aligner
    let aligner = IntegratedAligner::new(align_config)
        .with_plane_sweep(sweep_config);

    // Open output PAF writer
    let mut output = std::fs::File::create(output_file)?;
    use std::io::Write;

    println!("Processing queries with integrated FastGA + plane sweep...");

    // Process all queries, applying both FastGA plane sweep and SweepGA filters
    let total = aligner.process_all_queries(
        query_file,
        target_file,
        |alignment| {
            // SweepGA's per-alignment filter could go here
            // For example: check if alignment is in acceptable regions
            true
        },
        |query_name, alignments| {
            println!("Query {}: {} alignments after plane sweep", query_name, alignments.len());

            // Now apply SweepGA's sophisticated filtering
            // This is where you'd integrate with paf_filter.rs
            let filtered = apply_sweepga_filters(alignments);

            // Write filtered alignments to output
            for aln in filtered {
                writeln!(output, "{}", aln.to_paf()).ok();
            }
        },
    )?;

    println!("Processed {} alignments total", total);
    Ok(())
}

/// Placeholder for SweepGA's existing filtering logic
fn apply_sweepga_filters(alignments: Vec<fastga_rs::Alignment>) -> Vec<fastga_rs::Alignment> {
    // This is where SweepGA's existing logic would go:
    // - Scaffold creation
    // - Transitive chaining
    // - Rescue phase
    // - etc.
    alignments
}

/// Example showing memory-efficient streaming for huge genomes
fn sweepga_streaming_example(
    query_file: &Path,
    target_file: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use fastga_rs::integrated::AlignmentStream;

    let aligner = IntegratedAligner::new(Config::default())
        .with_plane_sweep(PlaneSweepConfig {
            max_per_query: 50,  // Very aggressive for huge genomes
            max_per_target: 1,
            min_identity: 0.95,
            min_length: 10000,
            max_overlap: 0.3,
        });

    // Create streaming iterator
    let stream = AlignmentStream::new(aligner, query_file, target_file)?;

    // Process alignments one by one without loading all into memory
    for alignment in stream {
        // Process each alignment immediately
        // No need to store everything in memory!
        if should_keep_alignment(&alignment) {
            println!("{}", alignment.to_paf());
        }
    }

    Ok(())
}

fn should_keep_alignment(alignment: &fastga_rs::Alignment) -> bool {
    // SweepGA's filtering logic
    alignment.identity() > 0.95
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Example usage
    let query_file = Path::new("data/queries.fasta");
    let target_file = Path::new("data/targets.fasta");
    let output_file = Path::new("filtered.paf");

    if query_file.exists() && target_file.exists() {
        sweepga_main_with_integrated_aligner(query_file, target_file, output_file)?;
    } else {
        println!("Example files not found. This shows how to integrate FastGA-rs into SweepGA.");
        println!("\nKey benefits:");
        println!("1. No intermediate PAF files for repetitive genomes");
        println!("2. Query-by-query processing ensures fair distribution");
        println!("3. Immediate plane sweep filtering reduces memory usage");
        println!("4. Direct integration - FastGA is embedded, no external binaries");
    }

    Ok(())
}