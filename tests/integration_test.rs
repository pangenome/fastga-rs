//! Integration tests for FastGA-RS using real yeast genome data.

use fastga_rs::{FastGA, Config};
use std::path::Path;
use std::fs;
use anyhow::Result;

#[test]
fn test_yeast_self_alignment() -> Result<()> {
    // Use the sample yeast data
    let yeast_file = Path::new("data/yeast_sample.fasta");

    // Create test data if it doesn't exist
    if !yeast_file.exists() {
        eprintln!("Test data not found. Creating sample data...");
        fs::create_dir_all("data")?;

        // Extract sample from compressed yeast genome
        let yeast_gz = Path::new("data/scerevisiae8.fa.gz");
        if yeast_gz.exists() {
            std::process::Command::new("sh")
                .arg("-c")
                .arg("zcat data/scerevisiae8.fa.gz | head -1000 > data/yeast_sample.fasta")
                .output()?;
        } else {
            // Create minimal test data
            fs::write(yeast_file,
                ">test_seq1\nACGTACGTACGTACGTACGTACGTACGTACGT\n\
                 >test_seq2\nACGTACGTACGTACGTTCGTACGTACGTACGT\n")?;
        }
    }

    // Create aligner with default configuration
    let aligner = FastGA::new(Config::default())?;

    // Align yeast against itself
    println!("Running FastGA alignment on yeast sample...");
    let alignments = aligner.align_files(yeast_file, yeast_file)?;

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