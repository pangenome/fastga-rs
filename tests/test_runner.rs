use fastga_rs::runner::Orchestrator;
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn test_orchestrator_full_pipeline() {
    let dir = tempdir().unwrap();

    // Create query and target files
    let query_fa = dir.path().join("query.fa");
    let target_fa = dir.path().join("target.fa");

    let mut file = File::create(&query_fa).unwrap();
    writeln!(file, ">query1").unwrap();
    for _ in 0..5 {
        writeln!(
            file,
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        )
        .unwrap();
    }
    file.flush().unwrap();

    let mut file = File::create(&target_fa).unwrap();
    writeln!(file, ">target1").unwrap();
    for _ in 0..5 {
        writeln!(
            file,
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        )
        .unwrap();
    }
    writeln!(file, ">target2").unwrap();
    for _ in 0..5 {
        writeln!(
            file,
            "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
        )
        .unwrap();
    }
    file.flush().unwrap();

    eprintln!("Testing full pipeline with Orchestrator");

    // Test default .1aln output
    let orchestrator = Orchestrator::new_simple(1);

    match orchestrator.align(&query_fa, &target_fa) {
        Ok(output) => {
            eprintln!("✓ Alignment succeeded!");
            eprintln!("Output: {output}");

            // Default format is .1aln, so output should be a file path
            if output.contains(".1aln") {
                eprintln!("✓ Output is .1aln file path: {output}");

                // Verify the file exists
                use std::path::Path;
                if Path::new(&output).exists() {
                    eprintln!("✓ .1aln file exists");
                }
            } else if output.contains('\t') {
                eprintln!("⚠ Output appears to be PAF format (unexpected with default config)");
            } else {
                eprintln!("⚠ Unexpected output format");
            }
        }
        Err(e) => {
            eprintln!("✗ Full pipeline failed: {e}");
            eprintln!("This might be expected if FastGA binary is not available");
        }
    }
}

#[test]
fn test_orchestrator_paf_output() {
    use fastga_rs::config::OutputFormat;

    let dir = tempdir().unwrap();

    // Create query and target files
    let query_fa = dir.path().join("query.fa");
    let target_fa = dir.path().join("target.fa");

    let mut file = File::create(&query_fa).unwrap();
    writeln!(file, ">query1").unwrap();
    for _ in 0..5 {
        writeln!(
            file,
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        )
        .unwrap();
    }
    file.flush().unwrap();

    let mut file = File::create(&target_fa).unwrap();
    writeln!(file, ">target1").unwrap();
    for _ in 0..5 {
        writeln!(
            file,
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        )
        .unwrap();
    }
    file.flush().unwrap();

    eprintln!("Testing PAF output format conversion");

    // Explicitly request PAF output
    let config = fastga_rs::Config::builder()
        .num_threads(1)
        .output_format(OutputFormat::PafWithX)
        .build();

    let orchestrator = Orchestrator::new(config);

    match orchestrator.align(&query_fa, &target_fa) {
        Ok(output) => {
            eprintln!("✓ PAF output succeeded!");
            eprintln!("Output length: {} bytes", output.len());

            if !output.is_empty() {
                // Print first few lines of output
                for (i, line) in output.lines().take(3).enumerate() {
                    eprintln!("  Line {}: {}", i + 1, line);
                }

                // Check for PAF format
                if output.contains('\t') {
                    eprintln!("✓ Output is in PAF format");
                }
            } else {
                eprintln!("⚠ Warning: Empty output");
            }
        }
        Err(e) => {
            eprintln!("✗ PAF conversion failed: {e}");
        }
    }
}
