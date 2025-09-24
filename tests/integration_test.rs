//! Integration tests for FastGA-RS using real yeast genome data.

use fastga_rs::{FastGA, Config};
use std::path::Path;
use std::fs;
use anyhow::Result;
use tempfile::TempDir;

#[test]
fn test_chrV_self_alignment() -> Result<()> {
    // Use chromosome V test data with all strains
    let chrV_file = Path::new("data/cerevisiae.chrV.fa.gz");

    if !chrV_file.exists() {
        eprintln!("Test data not found. Skipping test.");
        return Ok(());
    }

    // Decompress to temp directory to avoid GDB conflicts
    let temp_dir = tempfile::TempDir::new()?;
    let temp_chrV = temp_dir.path().join("chrV.fasta");

    // Decompress the file
    use std::process::Command;
    Command::new("gunzip")
        .arg("-c")
        .arg(chrV_file)
        .output()
        .and_then(|output| fs::write(&temp_chrV, output.stdout))?

    // Create aligner with default configuration
    let aligner = FastGA::new(Config::default())?;

    // Align chromosome against itself
    println!("Running FastGA alignment on chromosome V...");
    let alignments = aligner.align_files(&temp_chrV, &temp_chrV)?;

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
fn test_chrV_cross_strain_alignment() -> Result<()> {
    // Test alignment using full chrV file with multiple strains
    let chrV_file = Path::new("data/cerevisiae.chrV.fa.gz");

    if !chrV_file.exists() {
        eprintln!("Test data not found. Skipping test.");
        return Ok(());
    }

    // Decompress to temp directory
    let temp_dir = TempDir::new()?;
    let temp_chrV = temp_dir.path().join("chrV.fasta");

    use std::process::Command;
    Command::new("gunzip")
        .arg("-c")
        .arg(chrV_file)
        .output()
        .and_then(|output| fs::write(&temp_chrV, output.stdout))?;

    let aligner = FastGA::new(Config::default())?;
    // Self-alignment will find all strain-vs-strain alignments
    let alignments = aligner.align_files(&temp_chrV, &temp_chrV)?;

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