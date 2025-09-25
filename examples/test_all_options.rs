use fastga_rs::{FastGA, Config, config::OutputFormat};
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    println!("=== Testing All FastGA Configuration Options ===\n");

    // Prepare test data
    let chrv_gz = Path::new("data/cerevisiae.chrV.fa.gz");
    if !chrv_gz.exists() {
        eprintln!("chrV data not found");
        return;
    }

    let output = Command::new("gunzip")
        .arg("-c")
        .arg(chrv_gz)
        .output()
        .unwrap();
    std::fs::write("test.fa", &output.stdout).unwrap();

    // Test 1: Default config (should impose NO restrictions)
    println!("Test 1: Default config (no restrictions)");
    let config = Config::default();
    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();
    println!("  Default: {} alignments", result.len());
    let default_count = result.len();
    assert!(result.len() > 1000, "Default should not filter alignments");

    // Test 2: With length filter
    println!("\nTest 2: With min_alignment_length = 1000");
    let config = Config::builder()
        .min_alignment_length(1000)
        .build();
    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();
    println!("  With length filter: {} alignments", result.len());
    assert!(result.len() < default_count, "Should filter short alignments");

    // Test 3: With identity filter
    println!("\nTest 3: With min_identity = 0.95");
    let config = Config::builder()
        .min_identity(0.95)
        .build();
    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();
    println!("  With identity filter: {} alignments", result.len());

    // Test 4: With soft masking
    println!("\nTest 4: With soft_masking = true");
    let config = Config::builder()
        .soft_masking(true)
        .build();
    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();
    println!("  With soft masking: {} alignments", result.len());

    // Test 5: Different output formats
    println!("\nTest 5: Testing output formats");
    for (format, name) in &[
        (OutputFormat::PafWithX, "PAF with X"),
        (OutputFormat::PafWithM, "PAF with ="),
        (OutputFormat::PafShort, "PAF short CS"),
    ] {
        let config = Config::builder()
            .output_format(*format)
            .min_alignment_length(5000) // To get fewer alignments for inspection
            .build();
        let aligner = FastGA::new(config).unwrap();
        let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();

        if let Some(first) = result.alignments.first() {
            println!("  {} format - CIGAR sample: {}",
                     name,
                     &first.cigar.chars().take(50).collect::<String>());
        }
    }

    // Test 6: Advanced seeding parameters
    println!("\nTest 6: Advanced seeding parameters");
    let config = Config::builder()
        .adaptive_seed_cutoff(20)
        .min_chain_coverage(0.5)
        .chain_start_threshold(100)
        .build();
    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();
    println!("  With advanced seeding: {} alignments", result.len());

    // Test 7: Logging and temp directory
    println!("\nTest 7: Logging and temp directory");
    let log_path = PathBuf::from("test.log");
    let temp_dir = PathBuf::from("/tmp");
    let config = Config::builder()
        .verbose(true)
        .log_file(log_path.clone())
        .temp_dir(temp_dir)
        .build();
    let aligner = FastGA::new(config).unwrap();
    let _result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();
    if log_path.exists() {
        println!("  ✓ Log file created: {} bytes", std::fs::metadata(&log_path).unwrap().len());
        std::fs::remove_file(&log_path).unwrap();
    }

    // Test 8: Keep intermediates
    println!("\nTest 8: Keep intermediates");
    let config = Config::builder()
        .keep_intermediates(true)
        .build();
    let aligner = FastGA::new(config).unwrap();
    let _result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();

    let gdb_file = Path::new("test.1gdb");
    if gdb_file.exists() {
        println!("  ✓ Intermediate GDB file kept");
    }

    // Test 9: Thread configuration
    println!("\nTest 9: Thread configuration");
    for threads in &[1, 4, 8] {
        let config = Config::builder()
            .num_threads(*threads)
            .min_alignment_length(5000) // For speed
            .build();
        let aligner = FastGA::new(config).unwrap();
        let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();
        println!("  {} threads: {} alignments", threads, result.len());
    }

    // Test 10: Verify tags are preserved
    println!("\nTest 10: Tag preservation");
    let config = Config::default();
    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(Path::new("test.fa"), Path::new("test.fa")).unwrap();

    // Write to PAF and check tags
    result.write_paf("test_tags.paf").unwrap();
    let paf_content = std::fs::read_to_string("test_tags.paf").unwrap();
    let first_line = paf_content.lines().next().unwrap();
    let fields: Vec<_> = first_line.split('\t').collect();

    println!("  PAF has {} fields (12 required + tags)", fields.len());
    if fields.len() > 12 {
        println!("  Tags present: {}", &fields[12..].join(" "));
    }

    println!("\n✅ All configuration options tested successfully!");
}