use fastga_rs::{FastGA, Config};
use std::path::Path;
use std::process::Command;

fn main() {
    // First check if we have chrV data
    let chrv_gz = Path::new("data/cerevisiae.chrV.fa.gz");
    if !chrv_gz.exists() {
        eprintln!("No chrV data found");
        return;
    }

    // Extract it
    let output = Command::new("gunzip")
        .arg("-c")
        .arg(chrv_gz)
        .output()
        .unwrap();

    std::fs::write("chrV.fasta", &output.stdout).unwrap();
    println!("Extracted chrV: {} bytes", output.stdout.len());

    // Now align with our library
    let config = Config::builder()
        .num_threads(4)
        .min_alignment_length(1000)
        .min_identity(0.95)
        .build();

    let aligner = FastGA::new(config).unwrap();
    println!("\nRunning alignment with our fork/exec orchestrator...");
    let result = aligner.align_files(Path::new("chrV.fasta"), Path::new("chrV.fasta")).unwrap();

    println!("\n=== Alignment Results ===");
    println!("✅ Total alignments: {}", result.len());

    // Check a few alignments
    for (i, alignment) in result.iter().take(3).enumerate() {
        println!("\nAlignment {}:", i + 1);
        println!("  Query: {} [{}-{}]", alignment.query_name, alignment.query_start, alignment.query_end);
        println!("  Target: {} [{}-{}]", alignment.target_name, alignment.target_start, alignment.target_end);
        println!("  Identity: {:.2}%", alignment.identity() * 100.0);
        println!("  Alignment block: {} bp", alignment.query_end - alignment.query_start);
        let cigar = &alignment.cigar;
        if !cigar.is_empty() {
            // Check for extended CIGAR (= and X)
            let has_equals = cigar.contains('=');
            let has_x = cigar.contains('X');
            println!("  CIGAR uses extended format: {} (has '=': {}, has 'X': {})",
                     has_equals || has_x, has_equals, has_x);
            println!("  CIGAR preview: {}...", &cigar[..cigar.len().min(100)]);
        }
    }

    // Verify we're getting the expected alignment count for chrV self-alignment
    if result.len() > 1000 && result.len() < 2000 {
        println!("\n✅ VERIFICATION PASSED: Got expected ~1400 alignments for chrV self-alignment");
    } else {
        println!("\n⚠️  WARNING: Unexpected alignment count (expected ~1400)");
    }
}