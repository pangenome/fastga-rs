//! Integration tests for FastGA-RS using real yeast genome data.

use anyhow::Result;
use fastga_rs::{Config, FastGA};
use std::fs;
use std::path::Path;
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
        .and_then(|output| fs::write(&temp_chrV, output.stdout))?;

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
    assert!(
        fields.len() >= 12,
        "PAF line should have at least 12 fields"
    );

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

    assert!(
        !alignments.is_empty(),
        "Should find alignments between strains"
    );

    // Track self-alignments by sequence name to find the best/longest one
    let mut self_coverage: std::collections::HashMap<String, Vec<_>> =
        std::collections::HashMap::new();
    let mut cross_strain_alignments = 0;

    for aln in &alignments.alignments {
        if aln.query_name == aln.target_name {
            // Self-alignment
            self_coverage
                .entry(aln.query_name.clone())
                .or_insert(Vec::new())
                .push(aln);
        } else {
            cross_strain_alignments += 1;
            let identity = aln.identity();
            println!(
                "Cross-strain alignment: {} vs {} = {:.2}%",
                aln.query_name,
                aln.target_name,
                identity * 100.0
            );
            assert!(
                identity > 0.70,
                "Different strains should share some similarity"
            );
        }
    }

    // Check self-alignment coverage for each strain
    for (strain, alns) in &self_coverage {
        // Find the main self-alignment (should cover most/all of the sequence)
        let main_aln = alns
            .iter()
            .max_by_key(|a| a.query_end - a.query_start)
            .expect("Should have at least one self-alignment");

        let coverage =
            (main_aln.query_end - main_aln.query_start) as f64 / main_aln.query_len as f64;
        println!(
            "Self-alignment for {}: coverage={:.1}%, identity={:.1}%",
            strain,
            coverage * 100.0,
            main_aln.identity() * 100.0
        );

        // Main self-alignment should have near-complete coverage and perfect identity
        assert!(
            coverage > 0.95,
            "{} should have >95% self-alignment coverage",
            strain
        );
        assert!(
            main_aln.identity() > 0.99,
            "{} self-alignment should be >99% identity",
            strain
        );
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
    assert!(aligner
        .align_files(Path::new("nonexistent1.fa"), Path::new("nonexistent2.fa"))
        .is_err());
}
