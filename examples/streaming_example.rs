//! Example demonstrating streaming alignment processing with FastGA.
//!
//! This example shows various ways to use the streaming API to process
//! alignments efficiently without storing them all in memory.

use anyhow::Result;
use fastga_rs::streaming::{align_streaming_simple, BestHitFilter, StreamingAligner};
use fastga_rs::{Alignment, Config};
use std::collections::HashMap;
use std::path::Path;

fn main() -> Result<()> {
    println!("FastGA Streaming API Examples\n");

    // Get test data paths
    let genome1 = Path::new("data/yeast_sample.fasta");
    let genome2 = Path::new("data/yeast_sample.fasta");

    if !genome1.exists() || !genome2.exists() {
        eprintln!("Test data not found. Please ensure data/yeast_sample.fasta exists.");
        eprintln!("You can create it with: zcat data/scerevisiae8.fa.gz | head -1000 > data/yeast_sample.fasta");
        return Ok(());
    }

    // Example 1: Simple streaming with counting
    example_simple_streaming(genome1, genome2)?;

    // Example 2: Filtering alignments
    example_filtering(genome1, genome2)?;

    // Example 3: Best hit per query
    example_best_hits(genome1, genome2)?;

    // Example 4: Statistics aggregation
    example_statistics(genome1, genome2)?;

    // Example 5: Query vs all targets
    example_query_vs_all(genome2)?;

    // Example 6: Stream to file
    example_stream_to_file(genome1, genome2)?;

    Ok(())
}

/// Example 1: Simple streaming alignment with counting
fn example_simple_streaming(genome1: &Path, genome2: &Path) -> Result<()> {
    println!("=== Example 1: Simple Streaming ===\n");

    let mut alignment_count = 0;
    let mut total_bases = 0;

    let stats = align_streaming_simple(genome1, genome2, |alignment| {
        alignment_count += 1;
        total_bases += alignment.query_end - alignment.query_start;

        // Process first 5 alignments
        if alignment_count <= 5 {
            println!(
                "Alignment {}: {} -> {} ({:.1}% identity)",
                alignment_count,
                alignment.query_name,
                alignment.target_name,
                alignment.identity() * 100.0
            );
        }

        true // Keep all alignments
    })?;

    println!("\nProcessed {} alignments", stats.total_alignments);
    println!("Total aligned bases: {total_bases}");
    println!();

    Ok(())
}

/// Example 2: Filtering alignments by various criteria
fn example_filtering(genome1: &Path, genome2: &Path) -> Result<()> {
    println!("=== Example 2: Filtered Streaming ===\n");

    let mut aligner = StreamingAligner::new(Config::default());

    // Add multiple filters
    aligner
        .filter_min_identity(0.95)
        .filter_min_length(100)
        .filter(|aln| {
            // Custom filter: exclude self-alignments
            !(aln.query_name == aln.target_name && aln.query_start == aln.target_start)
        });

    let mut high_quality = Vec::new();

    let stats = aligner.align_files(genome1, genome2, |alignment| {
        high_quality.push(alignment.clone());
        true
    })?;

    println!("Total alignments: {}", stats.total_alignments);
    println!("After filtering: {}", high_quality.len());
    println!("Filtered out: {}", stats.filtered_alignments);

    // Show some high-quality alignments
    for (i, aln) in high_quality.iter().take(3).enumerate() {
        println!("\nHigh-quality alignment {}:", i + 1);
        println!("  {} -> {}", aln.query_name, aln.target_name);
        println!("  Identity: {:.2}%", aln.identity() * 100.0);
        println!("  Length: {} bp", aln.query_end - aln.query_start);
    }

    println!();
    Ok(())
}

/// Example 3: Keep only best hit per query
fn example_best_hits(genome1: &Path, genome2: &Path) -> Result<()> {
    println!("=== Example 3: Best Hit Per Query ===\n");

    let mut best_filter = BestHitFilter::new();

    align_streaming_simple(genome1, genome2, |alignment| {
        best_filter.process(alignment);
        true
    })?;

    let best_alignments = best_filter.into_alignments();

    println!("Found best hits for {} queries", best_alignments.len());

    for aln in best_alignments.iter().take(5) {
        println!(
            "Best hit for {}: {} ({:.2}% identity)",
            aln.query_name,
            aln.target_name,
            aln.identity() * 100.0
        );
    }

    println!();
    Ok(())
}

/// Example 4: Aggregate statistics during streaming
fn example_statistics(genome1: &Path, genome2: &Path) -> Result<()> {
    println!("=== Example 4: Statistics Aggregation ===\n");

    use std::sync::{Arc, Mutex};

    let identity_histogram = Arc::new(Mutex::new(vec![0; 11])); // 0-10%, 10-20%, ..., 90-100%
    let length_distribution = Arc::new(Mutex::new(HashMap::new()));
    let query_counts = Arc::new(Mutex::new(HashMap::new()));

    let mut aligner = StreamingAligner::new(Config::default());

    let hist_clone = Arc::clone(&identity_histogram);
    let len_clone = Arc::clone(&length_distribution);
    let query_clone = Arc::clone(&query_counts);

    aligner.aggregate(move |alignment: &Alignment| {
        // Update identity histogram
        let bucket = (alignment.identity() * 10.0) as usize;
        let mut hist = hist_clone.lock().unwrap();
        if bucket < hist.len() {
            hist[bucket] += 1;
        }

        // Track alignment lengths
        let length_bucket = ((alignment.query_end - alignment.query_start) / 100) * 100;
        let mut len_dist = len_clone.lock().unwrap();
        *len_dist.entry(length_bucket).or_insert(0) += 1;

        // Count alignments per query
        let mut counts = query_clone.lock().unwrap();
        *counts.entry(alignment.query_name.clone()).or_insert(0) += 1;
    });

    let stats = aligner.align_files(
        genome1,
        genome2,
        |_| true, // Keep all
    )?;

    println!("Processed {} alignments\n", stats.total_alignments);

    // Display identity distribution
    println!("Identity distribution:");
    let hist = identity_histogram.lock().unwrap();
    for (i, count) in hist.iter().enumerate() {
        if *count > 0 {
            let start = i * 10;
            let end = (i + 1) * 10;
            println!("  {start}%-{end}%: {count} alignments");
        }
    }

    // Display length distribution
    println!("\nLength distribution:");
    let len_dist = length_distribution.lock().unwrap();
    let mut lengths: Vec<_> = len_dist.iter().collect();
    lengths.sort_by_key(|&(k, _)| k);
    for (bucket, count) in lengths.iter().take(5) {
        println!("  {}-{} bp: {} alignments", bucket, *bucket + 100, count);
    }

    // Display query statistics
    println!("\nQueries with most alignments:");
    let counts = query_counts.lock().unwrap();
    let mut query_vec: Vec<_> = counts.iter().collect();
    query_vec.sort_by_key(|&(_, count)| std::cmp::Reverse(*count));
    for (query, count) in query_vec.iter().take(3) {
        println!("  {query}: {count} alignments");
    }

    println!();
    Ok(())
}

/// Example 5: Align single query against all targets
fn example_query_vs_all(target_file: &Path) -> Result<()> {
    println!("=== Example 5: Query vs All ===\n");

    // Create a test query sequence
    let query = b">test_query\n\
                  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n\
                  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

    let mut aligner = StreamingAligner::new(
        Config::builder()
            .min_alignment_length(50)
            .min_identity(0.8)
            .build(),
    );

    let alignments = aligner.align_query_vs_all(query, target_file, |alignment| {
        println!(
            "Found alignment to {}: {:.1}% identity",
            alignment.target_name,
            alignment.identity() * 100.0
        );
        alignment.identity() > 0.85 // Keep high-identity hits
    })?;

    println!(
        "\nQuery aligned to {} targets with >85% identity",
        alignments.len()
    );
    println!();

    Ok(())
}

/// Example 6: Stream alignments directly to PAF file
fn example_stream_to_file(genome1: &Path, genome2: &Path) -> Result<()> {
    println!("=== Example 6: Stream to PAF File ===\n");

    use std::fs::File;
    use std::io::BufWriter;

    let output_file = "alignments.paf";
    let file = File::create(output_file)?;
    let writer = BufWriter::new(file);

    let stats = fastga_rs::streaming::stream_to_paf(genome1, genome2, Config::default(), writer)?;

    println!(
        "Streamed {} alignments to {}",
        stats.kept_alignments, output_file
    );

    // Show file size
    let metadata = std::fs::metadata(output_file)?;
    println!("Output file size: {} bytes", metadata.len());

    // Clean up
    std::fs::remove_file(output_file)?;

    println!();
    Ok(())
}
