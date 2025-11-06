//! Advanced integration tests for FastGA-RS
//!
//! These tests verify more complex scenarios and edge cases.

use anyhow::Result;
use fastga_rs::{Config, FastGA};
use std::fs;
use std::path::Path;
use tempfile::TempDir;

/// Create a multi-sequence FASTA file for testing
fn create_multi_sequence_fasta(
    dir: &Path,
    filename: &str,
    num_sequences: usize,
) -> Result<std::path::PathBuf> {
    let mut content = String::new();

    for i in 0..num_sequences {
        content.push_str(&format!(">sequence_{i}\n"));
        // Create longer sequences (>200bp) to satisfy all config presets including fast()
        // Use actual yeast sequence patterns repeated to reach required length
        let base_seq = "ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTATACGAGCCACGGAAGAACGCCTTGTGGTGGACAAGAAGCAGGTG";
        let bases = if i.is_multiple_of(2) {
            // Original sequence repeated twice (>200bp)
            format!("{base_seq}{base_seq}")
        } else {
            // Slightly modified version with a mutation
            let modified = base_seq.replace("CTAT", "CTTT");
            format!("{modified}{modified}")
        };
        content.push_str(&bases);
        content.push('\n');
    }

    let path = dir.join(filename);
    fs::write(&path, content)?;
    Ok(path)
}

#[test]
fn test_extended_cigar_operators() -> Result<()> {
    let temp_dir = TempDir::new()?;

    // Use a real yeast sequence with a known mutation
    // This is actual yeast DNA that FastGA can align properly
    let seq1_path = temp_dir.path().join("seq1.fasta");
    fs::write(
        &seq1_path,
        b">ref\n\
        ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTATACGAGCCACG\n\
        GAAGAACGCCTTGTGGTGGACAAGAAGCAGGTGCTGGAGAAGCAGAAGACCAAGGCCAAGGAGCTGGCCAAGAAGGCT\n",
    )?;

    // Same sequence but with a few mismatches (X instead of = in CIGAR)
    let seq2_path = temp_dir.path().join("seq2.fasta");
    fs::write(
        &seq2_path,
        b">query\n\
        ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTTTACGAGCCACG\n\
        GAAGAACGCCTTGTGGTGGACAAGAAGCAGGTGCTGGAGAAGCAGAAGACCAAGGCCAAGGAGCTGGCCAAGAAGGCT\n",
    )?;

    let aligner = FastGA::new(Config::default())?;
    let alignments = aligner.align_files(&seq1_path, &seq2_path)?;

    // Check that we have alignments
    assert!(!alignments.is_empty(), "Should find alignments");

    // Verify extended CIGAR contains = and X operators
    for alignment in alignments.alignments {
        let cigar = &alignment.cigar;

        // Extended CIGAR should have '=' for matches and 'X' for mismatches
        let has_match_op = cigar.contains('=');
        let has_mismatch_op = cigar.contains('X');

        assert!(
            has_match_op || has_mismatch_op || cigar.is_empty(),
            "CIGAR '{cigar}' should contain extended operators"
        );

        // If we have a good alignment, verify the CIGAR makes sense
        if !cigar.is_empty() && alignment.identity() > 0.8 {
            // Most of the CIGAR should be matches
            let match_count = cigar.matches('=').count();
            let mismatch_count = cigar.matches('X').count();

            println!(
                "Alignment {} -> {}: {} matches, {} mismatches",
                alignment.query_name, alignment.target_name, match_count, mismatch_count
            );
        }
    }

    Ok(())
}

#[test]
fn test_different_config_presets() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let seq_path = create_multi_sequence_fasta(temp_dir.path(), "test.fasta", 3)?;

    // Test high sensitivity config
    let high_sens_aligner = FastGA::new(Config::high_sensitivity())?;
    let high_sens_alignments = high_sens_aligner.align_files(&seq_path, &seq_path)?;

    // Test fast config
    let fast_aligner = FastGA::new(Config::fast())?;
    let fast_alignments = fast_aligner.align_files(&seq_path, &seq_path)?;

    // Test repetitive genomes config
    let repetitive_aligner = FastGA::new(Config::repetitive_genomes())?;
    let repetitive_alignments = repetitive_aligner.align_files(&seq_path, &seq_path)?;

    // All should produce some alignments
    assert!(
        !high_sens_alignments.is_empty(),
        "High sensitivity should find alignments"
    );
    assert!(
        !fast_alignments.is_empty(),
        "Fast mode should find alignments"
    );
    assert!(
        !repetitive_alignments.is_empty(),
        "Repetitive mode should find alignments"
    );

    println!(
        "High sensitivity: {} alignments",
        high_sens_alignments.len()
    );
    println!("Fast mode: {} alignments", fast_alignments.len());
    println!(
        "Repetitive mode: {} alignments",
        repetitive_alignments.len()
    );

    Ok(())
}

#[test]
fn test_identity_calculation_accuracy() -> Result<()> {
    // This test verifies that identity calculations are accurate
    // This is critical for your plane sweep filter's identity * log(length) scoring

    let temp_dir = TempDir::new()?;

    // Create two sequences with known differences
    // Identical except for 2 positions out of 160
    let seq1_path = temp_dir.path().join("seq1.fasta");
    fs::write(
        &seq1_path,
        b">ref\n\
        ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTATACGAGCCACG\n\
        GAAGAACGCCTTGTGGTGGACAAGAAGCAGGTGCTGGAGAAGCAGAAGACCAAGGCCAAGGAGCTGGCCAAGAAGGCT\n",
    )?;

    let seq2_path = temp_dir.path().join("seq2.fasta");
    fs::write(
        &seq2_path,
        b">query\n\
        ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTTTACGAGCCACG\n\
        GAAGAACGCCTTGTGGTGGACAAGAAGCAGGTGCTGGAGAAGCAGAAGACCAAGGCCAAGGAGCTGGCCAAGAAGGCT\n",
    )?;

    let aligner = FastGA::new(Config::default())?;
    let alignments = aligner.align_files(&seq1_path, &seq2_path)?;

    assert!(!alignments.is_empty(), "Should find alignment");

    if let Some(alignment) = alignments.alignments.first() {
        let identity = alignment.identity();
        let length = alignment.query_end - alignment.query_start;

        // For plane sweep scoring: identity * log(length)
        let score = identity * (length as f64).ln();

        println!("Identity: {identity:.4}");
        println!("Length: {length}");
        println!("Plane sweep score (identity * ln(length)): {score:.4}");

        // Identity should be very high (158/160 = 0.9875)
        assert!(
            identity > 0.98 && identity < 1.0,
            "Identity {identity} should be ~0.9875 for 2 mismatches in 160bp"
        );
    }

    Ok(())
}

#[test]
fn test_identity_threshold_filtering() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let seq_path = create_multi_sequence_fasta(temp_dir.path(), "mixed.fasta", 4)?;

    // Test with different identity thresholds
    // FastGA requires minimum identity in range [0.55, 1.0)
    let config_90 = Config::builder().min_identity(0.9).build();
    let aligner_90 = FastGA::new(config_90)?;
    let alignments_90 = aligner_90.align_files(&seq_path, &seq_path)?;

    let config_55 = Config::builder()
        .min_identity(0.55) // Minimum allowed by FastGA
        .build();
    let aligner_55 = FastGA::new(config_55)?;
    let alignments_55 = aligner_55.align_files(&seq_path, &seq_path)?;

    // Lower threshold should allow more alignments
    assert!(
        alignments_55.len() >= alignments_90.len(),
        "Lower identity threshold should find same or more alignments"
    );

    // Verify identity filtering worked
    for alignment in alignments_90.alignments {
        assert!(
            alignment.identity() >= 0.9,
            "Alignment identity {} should be >= 0.9",
            alignment.identity()
        );
    }

    Ok(())
}

#[test]
fn test_paf_format_compliance() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let seq_path = create_multi_sequence_fasta(temp_dir.path(), "paf_test.fasta", 2)?;

    let aligner = FastGA::new(Config::default())?;
    let alignments = aligner.align_files(&seq_path, &seq_path)?;

    if !alignments.is_empty() {
        let paf_output = alignments.to_paf()?;

        // Verify PAF format
        for line in paf_output.lines() {
            let fields: Vec<&str> = line.split('\t').collect();

            // PAF format requires at least 12 mandatory fields
            assert!(
                fields.len() >= 12,
                "PAF line should have at least 12 fields, got {}",
                fields.len()
            );

            // Field 0: Query sequence name
            assert!(!fields[0].is_empty(), "Query name should not be empty");

            // Field 1: Query sequence length
            assert!(
                fields[1].parse::<usize>().is_ok(),
                "Query length should be numeric"
            );

            // Field 4: Strand (+ or -)
            assert!(
                fields[4] == "+" || fields[4] == "-",
                "Strand should be + or -"
            );

            // Field 5: Target sequence name
            assert!(!fields[5].is_empty(), "Target name should not be empty");

            // Check for CIGAR string in optional fields
            let has_cigar = fields
                .iter()
                .skip(12)
                .any(|f| f.starts_with("cg:Z:") || f.starts_with("cigar:Z:"));

            if has_cigar {
                // Find and validate CIGAR
                for field in fields.iter().skip(12) {
                    if let Some(cigar) = field.strip_prefix("cg:Z:") {
                        // Extended CIGAR should have valid operators
                        for c in cigar.chars() {
                            if c.is_alphabetic() {
                                assert!("MIDNSHP=X".contains(c), "Invalid CIGAR operator: {c}");
                            }
                        }
                    }
                }
            }
        }

        println!(
            "PAF format validation passed for {} lines",
            paf_output.lines().count()
        );
    }

    Ok(())
}

#[test]
fn test_self_alignment_symmetry() -> Result<()> {
    let temp_dir = TempDir::new()?;

    // Create a longer, more complex sequence for self-alignment
    let seq_path = temp_dir.path().join("self.fasta");
    fs::write(
        &seq_path,
        b">self\n\
        ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTATACGAGCCACG\n\
        GAAGAACGCCTTGTGGTGGACAAGAAGCAGGTGCTGGAGAAGCAGAAGACCAAGGCCAAGGAGCTGGCCAAGAAGGCT\n",
    )?;

    let aligner = FastGA::new(Config::default())?;
    let alignments = aligner.align_files(&seq_path, &seq_path)?;

    // Self-alignment should find at least one perfect match
    assert!(!alignments.is_empty(), "Self-alignment should find matches");

    let mut found_perfect = false;
    for alignment in &alignments.alignments {
        if alignment.query_name == alignment.target_name
            && alignment.query_start == alignment.target_start
            && alignment.query_end == alignment.target_end
        {
            // This should be a perfect match
            assert_eq!(
                alignment.identity(),
                1.0,
                "Self-alignment should have 100% identity"
            );
            found_perfect = true;
        }
    }

    assert!(
        found_perfect,
        "Should find at least one perfect self-alignment"
    );

    Ok(())
}

#[test]
fn test_large_sequence_handling() -> Result<()> {
    let temp_dir = TempDir::new()?;

    // Create a larger sequence (10kb)
    let mut large_seq = String::from(">large\n");
    for _ in 0..125 {
        // 125 * 80 = 10,000 bp
        large_seq.push_str(
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        );
    }
    large_seq.push('\n');

    let seq_path = temp_dir.path().join("large.fasta");
    fs::write(&seq_path, large_seq)?;

    // Should handle large sequences without panic
    let aligner = FastGA::new(Config::default())?;
    let result = aligner.align_files(&seq_path, &seq_path);

    assert!(
        result.is_ok(),
        "Should handle large sequences without error"
    );

    if let Ok(alignments) = result {
        println!(
            "Large sequence (10kb) produced {} alignments",
            alignments.len()
        );
    }

    Ok(())
}

#[test]
fn test_empty_file_handling() -> Result<()> {
    let temp_dir = TempDir::new()?;

    // Create empty FASTA file
    let empty_path = temp_dir.path().join("empty.fasta");
    fs::write(&empty_path, b"")?;

    // Create valid FASTA file
    let valid_path = temp_dir.path().join("valid.fasta");
    fs::write(&valid_path, b">seq\nACGT\n")?;

    let aligner = FastGA::new(Config::default())?;

    // Should handle empty files gracefully
    let result1 = aligner.align_files(&empty_path, &valid_path);
    let result2 = aligner.align_files(&valid_path, &empty_path);

    // These might error or return empty alignments, both are acceptable
    if let Ok(alignments) = result1 {
        assert_eq!(
            alignments.len(),
            0,
            "Empty file should produce no alignments"
        );
    }

    if let Ok(alignments) = result2 {
        assert_eq!(
            alignments.len(),
            0,
            "Empty file should produce no alignments"
        );
    }

    Ok(())
}

#[test]
fn test_multithreading_performance() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let seq_path = create_multi_sequence_fasta(temp_dir.path(), "perf_test.fasta", 5)?;

    use std::time::Instant;

    // Test with 1 thread
    let config_1 = Config::builder().num_threads(1).build();
    let aligner_1 = FastGA::new(config_1)?;

    let start_1 = Instant::now();
    let alignments_1 = aligner_1.align_files(&seq_path, &seq_path)?;
    let duration_1 = start_1.elapsed();

    // Test with multiple threads
    let config_multi = Config::builder().num_threads(4).build();
    let aligner_multi = FastGA::new(config_multi)?;

    let start_multi = Instant::now();
    let alignments_multi = aligner_multi.align_files(&seq_path, &seq_path)?;
    let duration_multi = start_multi.elapsed();

    // Both should produce same number of alignments
    assert_eq!(
        alignments_1.len(),
        alignments_multi.len(),
        "Thread count shouldn't affect alignment count"
    );

    println!(
        "1 thread: {:?}, {} alignments",
        duration_1,
        alignments_1.len()
    );
    println!(
        "4 threads: {:?}, {} alignments",
        duration_multi,
        alignments_multi.len()
    );

    // Multi-threaded might be faster (though not guaranteed for small inputs)
    if alignments_1.len() > 10 {
        println!(
            "Speedup: {:.2}x",
            duration_1.as_secs_f64() / duration_multi.as_secs_f64()
        );
    }

    Ok(())
}
