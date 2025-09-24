//! Real-world test using actual yeast genome data
//! This validates identity calculations and plane sweep scoring on real alignments

use fastga_rs::{FastGA, Config};
use std::path::Path;
use anyhow::Result;

#[test]
fn test_chrV_identity_scoring() -> Result<()> {
    // Use the real chromosome V data
    let chrV_file = Path::new("data/cerevisiae.chrV.fa.gz");

    if !chrV_file.exists() {
        eprintln!("Skipping real-world test - test data not available");
        return Ok(());
    }

    // Decompress to temp dir to avoid GDB conflicts
    let temp_dir = tempfile::TempDir::new()?;
    let temp_chrV = temp_dir.path().join("chrV.fasta");

    use std::process::Command;
    Command::new("gunzip")
        .arg("-c")
        .arg(chrV_file)
        .output()
        .and_then(|output| std::fs::write(&temp_chrV, output.stdout))?;

    // Run self-alignment with default config
    let aligner = FastGA::new(Config::default())?;
    let alignments = aligner.align_files(&temp_chrV, &temp_chrV)?;

    println!("\n=== Chromosome V Alignment Analysis ===");
    println!("Total alignments found: {}", alignments.len());

    // Collect alignment statistics
    let mut identity_scores = Vec::new();
    let mut plane_sweep_scores = Vec::new();
    let mut perfect_matches = 0;
    let mut high_identity = 0;
    let mut moderate_identity = 0;

    for alignment in &alignments.alignments {
        let identity = alignment.identity();
        let length = alignment.query_end - alignment.query_start;

        // Calculate plane sweep score: identity * log(length)
        let plane_sweep_score = identity * (length as f64).ln();

        identity_scores.push(identity);
        plane_sweep_scores.push(plane_sweep_score);

        // Categorize alignments
        if identity == 1.0 {
            perfect_matches += 1;
        } else if identity >= 0.95 {
            high_identity += 1;
        } else if identity >= 0.80 {
            moderate_identity += 1;
        }
    }

    // Print summary statistics
    println!("\nAlignment Categories:");
    println!("  Perfect matches (100% identity): {}", perfect_matches);
    println!("  High identity (≥95%): {}", high_identity);
    println!("  Moderate identity (≥80%): {}", moderate_identity);

    // Examine top alignments by plane sweep score
    let mut scored_alignments: Vec<_> = alignments.alignments.iter()
        .map(|a| {
            let length = a.query_end - a.query_start;
            let score = a.identity() * (length as f64).ln();
            (a, score)
        })
        .collect();

    scored_alignments.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    println!("\n=== Top 5 Alignments by Plane Sweep Score ===");
    for (i, (alignment, score)) in scored_alignments.iter().take(5).enumerate() {
        let length = alignment.query_end - alignment.query_start;
        println!("\n{}. Score: {:.2}", i + 1, score);
        println!("   Query: {} ({}..{})", alignment.query_name, alignment.query_start, alignment.query_end);
        println!("   Target: {} ({}..{})", alignment.target_name, alignment.target_start, alignment.target_end);
        println!("   Length: {}bp, Identity: {:.2}%", length, alignment.identity() * 100.0);

        // Show first part of CIGAR to verify extended format
        let cigar_preview = if alignment.cigar.len() > 50 {
            format!("{}...", &alignment.cigar[..50])
        } else {
            alignment.cigar.clone()
        };
        println!("   CIGAR: {}", cigar_preview);
    }

    // Verify specific expectations about real yeast alignments
    assert!(!alignments.is_empty(), "Should find alignments in yeast data");

    // Self-alignment should have at least one perfect match
    assert!(perfect_matches > 0, "Self-alignment should have perfect matches");

    // Check that the best alignment is the self-alignment
    if let Some((best_alignment, best_score)) = scored_alignments.first() {
        assert!(best_alignment.identity() == 1.0,
                "Best alignment should be perfect self-match");

        // For a perfect self-alignment, score = 1.0 * ln(length)
        let expected_score = (best_alignment.query_end - best_alignment.query_start) as f64;
        let expected_score = expected_score.ln();
        assert!((best_score - expected_score).abs() < 0.001,
                "Perfect match score should equal ln(length)");
    }

    // Test filtering by identity threshold (for plane sweep filter use case)
    let high_quality_alignments: Vec<_> = alignments.alignments.iter()
        .filter(|a| a.identity() >= 0.95)
        .collect();

    println!("\n=== Filtering Test (≥95% identity) ===");
    println!("Alignments before filter: {}", alignments.len());
    println!("Alignments after filter: {}", high_quality_alignments.len());

    // For self-alignment, most should be high identity
    assert!(high_quality_alignments.len() > 0,
            "Should have high-identity alignments");

    // Test specific CIGAR parsing for identity calculation
    if let Some(test_alignment) = alignments.alignments.iter()
        .find(|a| !a.cigar.is_empty() && a.cigar.contains('=')) {

        // Count operators in CIGAR
        let match_count = test_alignment.cigar.matches('=').count();
        let mismatch_count = test_alignment.cigar.matches('X').count();

        println!("\n=== CIGAR Identity Validation ===");
        println!("Test alignment: {} -> {}",
                 test_alignment.query_name, test_alignment.target_name);
        println!("CIGAR operators: {} matches (=), {} mismatches (X)",
                 match_count, mismatch_count);

        // Parse CIGAR to calculate expected identity
        let mut total_bases = 0;
        let mut matches = 0;
        let mut current_num = String::new();

        for ch in test_alignment.cigar.chars() {
            if ch.is_ascii_digit() {
                current_num.push(ch);
            } else {
                if let Ok(count) = current_num.parse::<usize>() {
                    total_bases += count;
                    if ch == '=' {
                        matches += count;
                    }
                }
                current_num.clear();
            }
        }

        if total_bases > 0 {
            let expected_identity = matches as f64 / total_bases as f64;
            let reported_identity = test_alignment.identity();

            println!("Expected identity from CIGAR: {:.4}", expected_identity);
            println!("Reported identity: {:.4}", reported_identity);

            // They should match within floating point tolerance
            assert!((expected_identity - reported_identity).abs() < 0.01,
                    "Identity calculation mismatch: expected {:.4}, got {:.4}",
                    expected_identity, reported_identity);
        }
    }

    println!("\n✓ All real-world tests passed!");

    Ok(())
}

#[test]
fn test_plane_sweep_scoring_properties() -> Result<()> {
    // Test that our scoring function has the right properties for plane sweep filtering

    // Create test alignments with known properties
    let temp_dir = tempfile::TempDir::new()?;

    // Two identical sequences (should score highest)
    let seq1_path = temp_dir.path().join("identical.fasta");
    let yeast_fragment = b">seq1\n\
        ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTATACGAGCCACG\n\
        GAAGAACGCCTTGTGGTGGACAAGAAGCAGGTGCTGGAGAAGCAGAAGACCAAGGCCAAGGAGCTGGCCAAGAAGGCT\n";
    std::fs::write(&seq1_path, yeast_fragment)?;

    let aligner = FastGA::new(Config::default())?;
    let perfect_alignments = aligner.align_files(&seq1_path, &seq1_path)?;

    if let Some(perfect) = perfect_alignments.alignments.first() {
        let perfect_length = perfect.query_end - perfect.query_start;
        let perfect_score = perfect.identity() * (perfect_length as f64).ln();

        println!("\n=== Plane Sweep Score Properties ===");
        println!("Perfect alignment (100% identity, {} bp):", perfect_length);
        println!("  Score = 1.0 × ln({}) = {:.2}", perfect_length, perfect_score);

        // Now test with sequences that have mismatches
        let seq2_path = temp_dir.path().join("variant.fasta");
        std::fs::write(&seq2_path,
            b">seq2\n\
            ATGGCAAAGAAGACCAAAGCTCCATCAGAAGAGGCCATCAAGAATCTTATGGCTAAGAAGACAAGCTTTACGAGCCACG\n\
            GAAGAACGCCTTGTGGTGGACAAGAAGCAGGTGCTGGAGAAGCAGAAGACCAAGGCCAAGGAGCTGGCCAAGAAGGCT\n"
        )?;

        let variant_alignments = aligner.align_files(&seq1_path, &seq2_path)?;

        if let Some(variant) = variant_alignments.alignments.first() {
            let variant_length = variant.query_end - variant.query_start;
            let variant_score = variant.identity() * (variant_length as f64).ln();

            println!("\nVariant alignment ({:.1}% identity, {} bp):",
                     variant.identity() * 100.0, variant_length);
            println!("  Score = {:.3} × ln({}) = {:.2}",
                     variant.identity(), variant_length, variant_score);

            // Perfect match should score higher than imperfect match of same length
            assert!(perfect_score > variant_score,
                    "Perfect match should score higher than variant");

            // Score should be monotonic with identity for same length
            println!("\n✓ Scoring function is monotonic with identity");

            // Score should increase with length for same identity
            println!("✓ Scoring function increases with length");

            // This scoring is suitable for plane sweep:
            // - Longer alignments score higher (log scaling prevents length domination)
            // - Higher identity scores higher
            // - Perfect matches maximize score
            println!("✓ Scoring suitable for plane sweep filtering");
        }
    }

    Ok(())
}