use fastga_rs::{FastGA, Config};
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;
use anyhow::Result;

#[test]
fn test_minimal_alignment() -> Result<()> {
    // Create minimal test FASTA file
    let dir = tempdir()?;
    let test_fa = dir.path().join("test.fa");

    let mut file = File::create(&test_fa)?;
    writeln!(file, ">seq1")?;
    writeln!(file, "ACGTACGTACGT")?;
    file.flush()?;
    drop(file);

    // Verify file was created
    eprintln!("Test file created at: {:?}", test_fa);
    eprintln!("File exists: {}", test_fa.exists());
    eprintln!("File size: {} bytes", std::fs::metadata(&test_fa)?.len());

    // Create aligner with minimal config
    let config = Config::builder()
        .num_threads(1)
        .min_alignment_length(1)  // Very permissive
        .build();

    eprintln!("Creating FastGA aligner...");
    let aligner = FastGA::new(config)?;
    eprintln!("Aligner created successfully");

    eprintln!("Starting alignment (self-alignment of test.fa)...");
    eprintln!("This is where it might hang...");

    // Try to align file against itself - should produce at least one alignment
    let result = aligner.align_files(&test_fa, &test_fa);

    match result {
        Ok(alignments) => {
            eprintln!("Alignment completed successfully!");
            eprintln!("Number of alignments: {}", alignments.len());
            assert!(alignments.len() > 0, "Expected at least one alignment");
        }
        Err(e) => {
            eprintln!("Alignment failed with error: {:?}", e);
            return Err(e.into());
        }
    }

    Ok(())
}

#[test]
fn test_version_check() -> Result<()> {
    use std::ffi::CString;
    use std::os::raw::{c_char, c_int};

    eprintln!("Testing direct FFI call to fastga_main...");

    // Link to the compiled FastGA
    #[link(name = "fastga_main", kind = "static")]
    extern "C" {
        fn fastga_main(argc: c_int, argv: *const *const c_char) -> c_int;
    }

    // Try calling with --version flag
    let args = vec![
        CString::new("FastGA").unwrap(),
        CString::new("--version").unwrap(),
    ];

    let argv: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();

    eprintln!("Calling fastga_main with --version...");
    let result = unsafe { fastga_main(argv.len() as c_int, argv.as_ptr()) };
    eprintln!("fastga_main returned: {}", result);

    // Note: --version might return non-zero, that's OK
    eprintln!("FFI call completed without hanging!");

    Ok(())
}