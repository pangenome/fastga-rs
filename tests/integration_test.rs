//! Integration tests for FastGA-RS using real yeast genome data.

use fastga_rs::{FastGA, Config};
use std::path::Path;
use std::fs;
use anyhow::Result;
use tempfile::TempDir;

#[test]
fn test_chr1_self_alignment() -> Result<()> {
    // Use chromosome I test data
    let chr1_ref = Path::new("data/chr1_ref.fasta");

    if !chr1_ref.exists() {
        eprintln!("Test data not found. Skipping test.");
        return Ok(());
    }

    // Copy to temp directory to avoid GDB conflicts
    let temp_dir = tempfile::TempDir::new()?;
    let temp_chr1 = temp_dir.path().join("chr1.fasta");
    fs::copy(chr1_ref, &temp_chr1)?;

    // Create aligner with default configuration
    let aligner = FastGA::new(Config::default())?;

    // Align chromosome against itself
    println!("Running FastGA alignment on chromosome I...");
    let alignments = aligner.align_files(&temp_chr1, &temp_chr1)?;

    // Basic validation
    assert!(!alignments.is_empty(), "Should find at least one alignment");
    println!("Found {} alignments", alignments.len());

    // Check first alignment
    if let Some(first) = alignments.alignments.first() {
        println!("First alignment:");
        println!("  Query: {} ({} bp)", first.query_name, first.query_len);
        println!("  Target: {} ({} bp)", first.target_name, first.target_len);
        println!("  Identity: {:.2}%", first.identity() * 100.0);
        println!("  CIGAR: {}", &first.cigar[..50.min(first.cigar.len())]);

        // Verify extended CIGAR format
        assert!(
            first.cigar.contains('=') || first.cigar.contains('X') || first.cigar.is_empty(),
            "CIGAR should contain extended operators (= or X)"
        );
    }

    // Convert to PAF and verify format
    let paf_output = alignments.to_paf()?;
    assert!(!paf_output.is_empty(), "PAF output should not be empty");

    // Check PAF format
    let first_line = paf_output.lines().next().unwrap();
    let fields: Vec<&str> = first_line.split('\t').collect();
    assert!(fields.len() >= 12, "PAF line should have at least 12 fields");

    println!("\nPAF output (first line):");
    println!("{}", first_line);

    Ok(())
}

#[test]
fn test_chr1_cross_strain_alignment() -> Result<()> {
    // Test alignment between two different yeast strains
    let chr1_ref = Path::new("data/chr1_ref.fasta");
    let chr1_s288c = Path::new("data/chr1_s288c.fasta");

    if !chr1_ref.exists() || !chr1_s288c.exists() {
        eprintln!("Test data not found. Skipping test.");
        return Ok(());
    }

    // Copy to temp directory
    let temp_dir = TempDir::new()?;
    let temp_ref = temp_dir.path().join("chr1_ref.fasta");
    let temp_s288c = temp_dir.path().join("chr1_s288c.fasta");
    fs::copy(chr1_ref, &temp_ref)?;
    fs::copy(chr1_s288c, &temp_s288c)?;

    let aligner = FastGA::new(Config::default())?;
    let alignments = aligner.align_files(&temp_ref, &temp_s288c)?;

    assert!(!alignments.is_empty(), "Should find alignments between strains");

    // Check identity - should be high but not perfect between strains
    if let Some(first) = alignments.alignments.first() {
        let identity = first.identity();
        println!("Cross-strain identity: {:.2}%", identity * 100.0);
        assert!(identity > 0.9, "Strains should be similar");
        assert!(identity < 1.0, "Strains should have some differences");
    }

    Ok(())
}

#[test]
fn test_config_presets() {
    // Test that preset configurations can be created
    let _high_sens = FastGA::new(Config::high_sensitivity()).unwrap();
    let _fast = FastGA::new(Config::fast()).unwrap();
    let _repetitive = FastGA::new(Config::repetitive_genomes()).unwrap();
}

#[test]
fn test_custom_config() {
    let config = Config::builder()
        .min_alignment_length(200)
        .min_identity(0.85)
        .num_threads(2)
        .build();

    let aligner = FastGA::new(config).unwrap();

    // This would test with actual alignment if we had test files
    // For now, just verify construction works
    assert!(aligner.align_files(
        Path::new("nonexistent1.fa"),
        Path::new("nonexistent2.fa")
    ).is_err());
}