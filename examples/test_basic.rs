use fastga_rs::{FastGA, Config};
use std::fs;
use std::path::Path;

fn main() -> anyhow::Result<()> {
    // Create test sequences
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".repeat(10);

    fs::write("test1.fa", format!(">seq1\n{}\n", seq))?;
    fs::write("test2.fa", format!(">seq2\n{}\n", seq))?;

    println!("Created test files");

    let config = Config::builder()
        .num_threads(1)
        .min_alignment_length(50)
        .build();

    println!("Running alignment...");
    let aligner = FastGA::new(config)?;
    let alignments = aligner.align_files(Path::new("test1.fa"), Path::new("test2.fa"))?;

    println!("Got {} alignments", alignments.len());

    if alignments.is_empty() {
        println!("ERROR: No alignments produced!");
    } else {
        for a in alignments.iter().take(5) {
            println!("  {} -> {}: {} bp, {:.1}% identity",
                a.query_name, a.target_name,
                a.query_end - a.query_start,
                a.identity() * 100.0);
        }
    }

    Ok(())
}