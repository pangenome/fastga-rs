//! Comprehensive tests for the streaming API.

use anyhow::Result;
use fastga_rs::streaming::{align_streaming_simple, BestHitFilter, StreamingAligner};
use fastga_rs::{Alignment, Config, FastGA};
use std::fs;
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Create test FASTA files using chromosome test data
fn create_test_sequences() -> Result<(tempfile::TempDir, std::path::PathBuf, std::path::PathBuf)> {
    let temp_dir = tempfile::TempDir::new()?;

    // Use chromosome V test data
    let chrV_file = Path::new("data/cerevisiae.chrV.fa.gz");

    if chrV_file.exists() {
        // Decompress to temp directory to avoid GDB conflicts
        let seq1_path = temp_dir.path().join("seq1.fasta");
        let seq2_path = temp_dir.path().join("seq2.fasta");

        // Decompress the file
        use std::process::Command;
        let output = Command::new("gunzip").arg("-c").arg(chrV_file).output()?;

        fs::write(&seq1_path, &output.stdout)?;
        fs::write(&seq2_path, &output.stdout)?;

        return Ok((temp_dir, seq1_path, seq2_path));
    }

    // Fallback if test data missing
    eprintln!("Warning: Test data not found, skipping test");
    let seq1_path = temp_dir.path().join("seq1.fasta");
    let seq2_path = temp_dir.path().join("seq2.fasta");

    // Create empty files to allow test to skip gracefully
    fs::write(&seq1_path, "")?;
    fs::write(&seq2_path, "")?;

    Ok((temp_dir, seq1_path, seq2_path))
}

#[test]
fn test_streaming_basic() -> Result<()> {
    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let mut alignments_collected = Vec::new();

    let stats = align_streaming_simple(&seq1, &seq2, |alignment| {
        alignments_collected.push(alignment);
        true // Keep all alignments
    })?;

    // We should get some alignments
    assert!(stats.total_alignments > 0);
    assert_eq!(stats.kept_alignments, alignments_collected.len());

    Ok(())
}

#[test]
fn test_streaming_with_filters() -> Result<()> {
    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let mut aligner = StreamingAligner::new(Config::default());

    // Add filters
    aligner
        .filter_min_identity(0.95) // High identity filter
        .filter_min_length(50); // Minimum length

    let mut high_quality_alignments = Vec::new();

    let stats = aligner.align_files(&seq1, &seq2, |alignment| {
        high_quality_alignments.push(alignment);
        true
    })?;

    // Check that filtering worked
    assert!(stats.filtered_alignments > 0 || stats.kept_alignments > 0);

    // All collected alignments should meet our criteria
    for aln in &high_quality_alignments {
        assert!(aln.identity() >= 0.95);
        assert!(aln.query_end - aln.query_start >= 50);
    }

    Ok(())
}

#[test]
fn test_query_filter() -> Result<()> {
    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let mut aligner = StreamingAligner::new(Config::default());

    // Filter only chr1 alignments
    aligner.filter_query(|name| name.contains("chr1"));

    let mut chr1_alignments = Vec::new();

    aligner.align_files(&seq1, &seq2, |alignment| {
        chr1_alignments.push(alignment.clone());
        true
    })?;

    // All alignments should be from chr1
    for aln in &chr1_alignments {
        assert!(aln.query_name.contains("chr1"));
    }

    Ok(())
}

#[test]
fn test_best_hit_filter() -> Result<()> {
    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let best_hits = Arc::new(Mutex::new(BestHitFilter::new()));
    let best_hits_clone = Arc::clone(&best_hits);

    let mut aligner = StreamingAligner::new(Config::default());

    aligner.align_files(&seq1, &seq2, move |alignment| {
        best_hits_clone.lock().unwrap().process(alignment);
        true
    })?;

    let best = Arc::try_unwrap(best_hits)
        .unwrap()
        .into_inner()
        .unwrap()
        .into_alignments();

    // Should have at most one alignment per query
    let unique_queries: std::collections::HashSet<_> = best.iter().map(|a| &a.query_name).collect();
    assert_eq!(unique_queries.len(), best.len());

    Ok(())
}

#[test]
fn test_streaming_callback_control() -> Result<()> {
    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let mut count = 0;
    let max_alignments = 2;

    let stats = align_streaming_simple(&seq1, &seq2, |_alignment| {
        count += 1;
        count <= max_alignments // Stop after max_alignments
    })?;

    // We should have processed some but kept only max_alignments
    assert!(stats.kept_alignments <= max_alignments);

    Ok(())
}

#[test]
fn test_aggregator() -> Result<()> {
    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let total_bases = Arc::new(Mutex::new(0usize));
    let total_bases_clone = Arc::clone(&total_bases);

    let mut aligner = StreamingAligner::new(Config::default());

    // Add an aggregator that counts total aligned bases
    aligner.aggregate(move |alignment: &Alignment| {
        let bases = alignment.query_end - alignment.query_start;
        *total_bases_clone.lock().unwrap() += bases;
    });

    aligner.align_files(&seq1, &seq2, |_| true)?;

    let total = *total_bases.lock().unwrap();
    assert!(total > 0, "Should have counted some aligned bases");

    Ok(())
}

#[test]
fn test_query_vs_all() -> Result<()> {
    let (_temp_dir, _seq1, seq2) = create_test_sequences()?;

    // Create a single query sequence
    let query = b">query\n\
                  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

    let mut aligner = StreamingAligner::new(Config::default());

    let alignments = aligner.align_query_vs_all(query, &seq2, |alignment| {
        // Only keep high-identity alignments
        alignment.identity() > 0.9
    })?;

    // All returned alignments should have the query name
    for aln in &alignments {
        assert_eq!(aln.query_name, "query");
        assert!(aln.identity() > 0.9);
    }

    Ok(())
}

#[test]
fn test_streaming_with_yeast_data() -> Result<()> {
    let yeast_file = Path::new("data/yeast_sample.fasta");

    if !yeast_file.exists() {
        eprintln!("Skipping yeast test - test data not available");
        return Ok(());
    }

    // Copy to temp directory to avoid GDB conflicts
    let temp_dir = tempfile::TempDir::new()?;
    let temp_yeast = temp_dir.path().join("yeast.fasta");
    fs::copy(yeast_file, &temp_yeast)?;

    let mut alignment_count = 0;
    let mut total_identity = 0.0;

    let stats = align_streaming_simple(
        &temp_yeast,
        &temp_yeast, // Self-alignment
        |alignment| {
            alignment_count += 1;
            total_identity += alignment.identity();

            // For self-alignment, same sequences should have 100% identity
            if alignment.query_name == alignment.target_name
                && alignment.query_start == alignment.target_start
            {
                assert_eq!(alignment.identity(), 1.0);
            }

            true
        },
    )?;

    assert!(alignment_count > 0);
    assert_eq!(alignment_count, stats.kept_alignments);

    let avg_identity = total_identity / alignment_count as f64;
    println!("Average identity: {:.2}%", avg_identity * 100.0);

    Ok(())
}

#[test]
fn test_paf_streaming() -> Result<()> {
    use std::io::Cursor;

    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let mut output = Cursor::new(Vec::new());

    let stats = fastga_rs::streaming::stream_to_paf(&seq1, &seq2, Config::default(), &mut output)?;

    let paf_output = String::from_utf8(output.into_inner())?;

    // Should have PAF format output
    assert!(!paf_output.is_empty());

    // Count lines to verify we got the right number of alignments
    let line_count = paf_output.lines().count();
    assert_eq!(line_count, stats.kept_alignments);

    // Verify PAF format (12+ tab-separated fields)
    for line in paf_output.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert!(
            fields.len() >= 12,
            "PAF line should have at least 12 fields"
        );
    }

    Ok(())
}

#[test]
fn test_multiple_filters_composition() -> Result<()> {
    let (_temp_dir, seq1, seq2) = create_test_sequences()?;

    let mut aligner = StreamingAligner::new(Config::default());

    // Add multiple filters that compose
    aligner
        .filter_min_identity(0.8)
        .filter_min_length(30)
        .filter_query(|name| !name.contains("chr2"))
        .filter_target(|name| name.contains("variant"))
        .filter(|aln| aln.mapping_quality >= 200);

    let mut filtered = Vec::new();

    aligner.align_files(&seq1, &seq2, |alignment| {
        filtered.push(alignment.clone());
        true
    })?;

    // Verify all filters were applied
    for aln in &filtered {
        assert!(aln.identity() >= 0.8);
        assert!(aln.query_end - aln.query_start >= 30);
        assert!(!aln.query_name.contains("chr2"));
        assert!(aln.target_name.contains("variant"));
        assert!(aln.mapping_quality >= 200);
    }

    Ok(())
}

#[test]
fn test_error_handling_invalid_file() {
    let mut aligner = StreamingAligner::new(Config::default());

    let result = aligner.align_files("nonexistent1.fa", "nonexistent2.fa", |_| true);

    assert!(result.is_err());
}
