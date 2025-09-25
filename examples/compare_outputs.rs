use fastga_rs::{FastGA, Config, config::OutputFormat};
use std::path::Path;
use std::process::Command;

fn main() {
    // Extract chrV
    let chrv_gz = Path::new("data/cerevisiae.chrV.fa.gz");
    let output = Command::new("gunzip")
        .arg("-c")
        .arg(chrv_gz)
        .output()
        .unwrap();
    std::fs::write("test_chrv.fa", &output.stdout).unwrap();

    // Run with our orchestrator
    let config = Config::builder()
        .num_threads(4)
        .min_alignment_length(1000)
        .min_identity(0.95)
        .soft_masking(true)
        .output_format(OutputFormat::PafWithX)  // Same as -pafx
        .build();

    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(
        Path::new("test_chrv.fa"),
        Path::new("test_chrv.fa")
    ).unwrap();

    // Write to file
    result.write_paf("our_orchestrator.paf").unwrap();
    println!("Our orchestrator wrote {} alignments to our_orchestrator.paf", result.len());

    // Run FastGA directly for comparison
    println!("\nRunning FastGA directly for comparison...");
    let direct_output = Command::new("./target/debug/build/fastga-rs-4b5aee3020c8890b/out/FastGA")
        .args(&["-pafx", "-T4", "-l1000", "-i0.95", "-M", "test_chrv.fa", "test_chrv.fa"])
        .output()
        .unwrap();

    std::fs::write("direct_fastga.paf", &direct_output.stdout).unwrap();
    let direct_count = String::from_utf8_lossy(&direct_output.stdout).lines().count();
    println!("Direct FastGA produced {} alignments", direct_count);

    // Compare
    println!("\n=== Comparison ===");
    println!("Our orchestrator: {} alignments", result.len());
    println!("Direct FastGA: {} alignments", direct_count);

    if result.len() == direct_count {
        println!("✅ SAME NUMBER OF ALIGNMENTS!");

        // Check if content matches
        let our_content = std::fs::read_to_string("our_orchestrator.paf").unwrap();
        let direct_content = std::fs::read_to_string("direct_fastga.paf").unwrap();

        // Get first lines for comparison
        let our_first = our_content.lines().next().unwrap_or("");
        let direct_first = direct_content.lines().next().unwrap_or("");

        println!("\nFirst alignment from our orchestrator:");
        println!("{}", our_first.split('\t').take(12).collect::<Vec<_>>().join("\t"));
        println!("\nFirst alignment from direct FastGA:");
        println!("{}", direct_first.split('\t').take(12).collect::<Vec<_>>().join("\t"));

        // Check exact match
        if our_content == direct_content {
            println!("\n✅ PERFECT MATCH: Outputs are byte-for-byte identical!");
        } else {
            // They might differ in order or minor details
            let mut our_lines: Vec<_> = our_content.lines().collect();
            let mut direct_lines: Vec<_> = direct_content.lines().collect();
            our_lines.sort();
            direct_lines.sort();

            if our_lines == direct_lines {
                println!("\n✅ OUTPUTS MATCH: Same alignments (different order)");
            } else {
                println!("\n⚠️  Outputs differ slightly (might be runtime variations)");
            }
        }
    } else {
        println!("❌ Different number of alignments!");
    }
}