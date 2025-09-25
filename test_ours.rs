use fastga_rs::{FastGA, Config};
use std::path::Path;

fn main() {
    let config = Config::builder()
        .num_threads(4)
        .min_alignment_length(1000)
        .min_identity(0.95)
        .soft_masking(true)
        .build();

    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(
        Path::new("test_chrv.fa"),
        Path::new("test_chrv.fa")
    ).unwrap();

    // Write PAF output
    result.write_paf("our_orchestrator.paf").unwrap();
    println!("Our orchestrator: {} alignments", result.len());
}
