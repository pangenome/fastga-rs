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
    let result = aligner.align_files(Path::new("chrV.fasta"), Path::new("chrV.fasta")).unwrap();
    
    println!("\n=== Alignment Results ===");
    println!("Total alignments: {}", result.len());
    
    // Check a few alignments
    for (i, alignment) in result.iter().take(3).enumerate() {
        println!("\nAlignment {}:", i + 1);
        println!("  Query: {} [{}-{}]", alignment.query_name, alignment.query_start, alignment.query_end);
        println!("  Target: {} [{}-{}]", alignment.target_name, alignment.target_start, alignment.target_end);
        println!("  Identity: {:.2}%", alignment.identity() * 100.0);
        println!("  Length: {} bp", alignment.alignment_length);
        if let Some(cigar) = &alignment.cigar {
            println!("  CIGAR preview: {}...", &cigar[..cigar.len().min(50)]);
        }
    }
}
