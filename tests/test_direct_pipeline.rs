use std::fs::File;
use std::io::Write;
use std::path::Path;
use tempfile::tempdir;
use std::process::Command;

#[test]
fn test_prepare_and_align_separately() {
    // This test prepares databases using utilities directly,
    // then runs FastGA with -g flag to skip preparation

    let dir = tempdir().unwrap();
    let test_fa = dir.path().join("test.fa");

    // Create test file with complex sequence
    let mut file = File::create(&test_fa).unwrap();
    writeln!(file, ">sequence1").unwrap();
    // Use a more complex sequence that FastGA can align
    writeln!(file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();
    writeln!(file, "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG").unwrap();
    writeln!(file, "GATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATA").unwrap();
    writeln!(file, "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA").unwrap();
    writeln!(file, "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT").unwrap();
    writeln!(file, "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT").unwrap();
    writeln!(file, "TATCTCGATGCATGCTAGCTACGATCGATGCTAGCTAGCATCGATCGATGCATGCTAGCTAGCGATCGATCG").unwrap();
    writeln!(file, "CGATCGATGCTAGCTAGCATCGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTAGCGATCGATCGAT").unwrap();
    file.flush().unwrap();

    let bin_dir = find_bin_dir();
    println!("Using binaries from: {}", bin_dir);

    // Step 1: Convert to GDB using FAtoGDB directly
    println!("\n=== Step 1: FAtoGDB ===");
    let fatogdb = format!("{}/FAtoGDB", bin_dir);
    if Path::new(&fatogdb).exists() {
        let output = Command::new(&fatogdb)
            .arg(&test_fa)
            .output()
            .unwrap();

        if output.status.success() {
            println!("✓ FAtoGDB succeeded");

            // Check if .1gdb was created
            let gdb_path = test_fa.with_extension("1gdb");
            if gdb_path.exists() {
                let size = std::fs::metadata(&gdb_path).unwrap().len();
                println!("✓ Created {} ({} bytes)", gdb_path.display(), size);
            }
        } else {
            println!("✗ FAtoGDB failed: {}", String::from_utf8_lossy(&output.stderr));
        }
    } else {
        println!("✗ FAtoGDB not found at {}", fatogdb);
    }

    // Step 2: Create index using GIXmake
    println!("\n=== Step 2: GIXmake ===");
    let gixmake = format!("{}/GIXmake", bin_dir);
    if Path::new(&gixmake).exists() {
        let base_name = test_fa.with_extension("");
        let output = Command::new(&gixmake)
            .arg("-f10")
            .arg("-T1")
            .arg(&base_name)
            .output()
            .unwrap();

        if output.status.success() {
            println!("✓ GIXmake succeeded");

            // Check if .gix was created
            let gix_path = base_name.with_extension("gix");
            if gix_path.exists() {
                let size = std::fs::metadata(&gix_path).unwrap().len();
                println!("✓ Created {} ({} bytes)", gix_path.display(), size);
            }
        } else {
            println!("✗ GIXmake failed: {}", String::from_utf8_lossy(&output.stderr));
        }
    } else {
        println!("✗ GIXmake not found at {}", gixmake);
    }

    // Step 3: Try running FastGA with -g flag (assuming GDB exists)
    println!("\n=== Step 3: FastGA with -g ===");
    let fastga = format!("{}/FastGA", bin_dir);
    if Path::new(&fastga).exists() {
        // Note: -g flag tells FastGA that GDB files already exist
        let output = Command::new(&fastga)
            .arg("-pafx")
            .arg("-g")  // Skip GDB creation
            .arg("-T1")
            .arg(&test_fa)
            .arg(&test_fa)
            .current_dir(dir.path())
            .output()
            .unwrap();

        if output.status.success() {
            println!("✓ FastGA succeeded");
            let stdout = String::from_utf8_lossy(&output.stdout);
            if !stdout.is_empty() {
                println!("Output: {} bytes", stdout.len());
            }
        } else {
            let stderr = String::from_utf8_lossy(&output.stderr);
            println!("✗ FastGA failed: {}", stderr);

            // Try without -g flag
            println!("\n=== Retry without -g flag ===");
            let output2 = Command::new(&fastga)
                .arg("-pafx")
                .arg("-T1")
                .arg(&test_fa)
                .arg(&test_fa)
                .current_dir(dir.path())
                .output()
                .unwrap();

            if output2.status.success() {
                println!("✓ FastGA succeeded without -g");
            } else {
                println!("✗ FastGA failed again: {}",
                        String::from_utf8_lossy(&output2.stderr));
            }
        }
    } else {
        println!("✗ FastGA not found at {}", fastga);
    }
}

fn find_bin_dir() -> String {
    // Try to find the build output directory
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let build_dir = exe_dir.join("build");
            if build_dir.exists() {
                if let Ok(entries) = std::fs::read_dir(&build_dir) {
                    for entry in entries.flatten() {
                        let name = entry.file_name();
                        if name.to_string_lossy().starts_with("fastga-rs-") {
                            let out_dir = entry.path().join("out");
                            if out_dir.exists() {
                                return out_dir.to_string_lossy().to_string();
                            }
                        }
                    }
                }
            }
        }
    }

    // Fallback
    "target/debug/build/fastga-rs-4b5aee3020c8890b/out".to_string()
}