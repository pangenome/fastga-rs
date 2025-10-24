//! Simple runner that uses the FastGA binary directly via subprocess
//! This avoids the FFI complexity and memory issues

use crate::error::{FastGAError, Result};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Find the FastGA binary
fn find_fastga_binary() -> Result<PathBuf> {
    // First try to find in the build output directories
    // This is where our compiled FastGA should be

    // Try to find in current exe's directory structure first
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            // Check relative to test binary location
            let search_dirs = vec![
                exe_dir.join("build"),
                exe_dir.parent().unwrap_or(exe_dir).join("build"),
            ];

            for base_dir in search_dirs {
                if base_dir.exists() {
                    for entry in std::fs::read_dir(&base_dir)
                        .unwrap_or_else(|_| std::fs::read_dir(".").unwrap())
                        .flatten()
                    {
                        let name = entry.file_name();
                        let name_str = name.to_string_lossy();
                        if name_str.starts_with("fastga-rs-") {
                            let fastga = entry.path().join("out/FastGA");
                            if fastga.exists() {
                                eprintln!("[FastGA] Using binary from build: {fastga:?}");
                                return Ok(fastga);
                            }
                        }
                    }
                }
            }
        }
    }

    // Check in OUT_DIR (from build)
    if let Ok(out_dir) = std::env::var("OUT_DIR") {
        let path = PathBuf::from(&out_dir).join("FastGA");
        if path.exists() {
            eprintln!("[FastGA] Using binary from OUT_DIR: {path:?}");
            return Ok(path);
        }
    }

    // Check in target/debug/build/*/out/
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    let base = PathBuf::from(&manifest_dir);

    // Try debug build
    let pattern = base.join("target/debug/build");
    if pattern.exists() {
        if let Ok(entries) = std::fs::read_dir(&pattern) {
            for entry in entries.flatten() {
                let name = entry.file_name();
                let name_str = name.to_string_lossy();
                if name_str.starts_with("fastga-rs-") {
                    let fastga = entry.path().join("out/FastGA");
                    if fastga.exists() {
                        eprintln!("[FastGA] Using binary from debug build: {fastga:?}");
                        return Ok(fastga);
                    }
                }
            }
        }
    }

    // Try release build
    let pattern = base.join("target/release/build");
    if pattern.exists() {
        if let Ok(entries) = std::fs::read_dir(&pattern) {
            for entry in entries.flatten() {
                let name = entry.file_name();
                let name_str = name.to_string_lossy();
                if name_str.starts_with("fastga-rs-") {
                    let fastga = entry.path().join("out/FastGA");
                    if fastga.exists() {
                        eprintln!("[FastGA] Using binary from release build: {fastga:?}");
                        return Ok(fastga);
                    }
                }
            }
        }
    }

    Err(FastGAError::Other(
        "FastGA binary not found. Please run 'cargo build' first.".to_string(),
    ))
}

/// Simple runner that directly executes FastGA
pub fn run_fastga_simple(
    query_path: &Path,
    target_path: &Path,
    num_threads: usize,
    min_length: usize,
    min_identity: Option<f64>,
) -> Result<String> {
    let fastga_path = find_fastga_binary()?;

    // Make sure the FastGA binary has all utilities in the same directory
    let bin_dir = fastga_path
        .parent()
        .ok_or_else(|| FastGAError::Other("Invalid FastGA path".to_string()))?;

    eprintln!("[FastGA] Binary directory: {bin_dir:?}");
    eprintln!("[FastGA] Checking for utilities...");

    // Check for required utilities
    for utility in &["FAtoGDB", "GIXmake"] {
        let util_path = bin_dir.join(utility);
        if !util_path.exists() {
            return Err(FastGAError::Other(format!(
                "Required utility {utility} not found in {bin_dir:?}"
            )));
        }
        eprintln!("[FastGA] Found utility: {util_path:?}");
    }

    // Build command
    let mut cmd = Command::new(&fastga_path);

    // Add arguments
    cmd.arg("-pafx"); // PAF with extended CIGAR
    cmd.arg(format!("-T{num_threads}"));
    cmd.arg(format!("-l{min_length}"));

    if let Some(identity) = min_identity {
        cmd.arg(format!("-i{identity:.2}"));
    }

    cmd.arg(query_path);
    cmd.arg(target_path);

    // Set ISOLATED PATH so FastGA can ONLY find its own utilities
    // This prevents FastGA from accidentally using system binaries
    cmd.env("PATH", bin_dir);

    eprintln!("[FastGA] Running command: {cmd:?}");

    // Execute and capture output
    let output = cmd
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| FastGAError::Other(format!("Failed to execute FastGA: {e}")))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        eprintln!("[FastGA] Command failed with stderr: {stderr}");
        return Err(FastGAError::FastGAExecutionFailed(stderr.to_string()));
    }

    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

/// Run FastGA with streaming output
pub fn run_fastga_streaming<F>(
    query_path: &Path,
    target_path: &Path,
    num_threads: usize,
    min_length: usize,
    min_identity: Option<f64>,
    mut callback: F,
) -> Result<()>
where
    F: FnMut(&str) -> bool,
{
    let fastga_path = find_fastga_binary()?;
    let bin_dir = fastga_path
        .parent()
        .ok_or_else(|| FastGAError::Other("Invalid FastGA path".to_string()))?;

    let mut cmd = Command::new(&fastga_path);

    cmd.arg("-pafx");
    cmd.arg(format!("-T{num_threads}"));
    cmd.arg(format!("-l{min_length}"));

    if let Some(identity) = min_identity {
        cmd.arg(format!("-i{identity:.2}"));
    }

    cmd.arg(query_path);
    cmd.arg(target_path);

    // Set ISOLATED PATH so FastGA can ONLY find its own utilities
    cmd.env("PATH", bin_dir);

    let mut child = cmd
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| FastGAError::Other(format!("Failed to spawn FastGA: {e}")))?;

    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| FastGAError::Other("Failed to capture stdout".to_string()))?;

    let reader = BufReader::new(stdout);

    for line in reader.lines() {
        let line = line.map_err(FastGAError::IoError)?;
        if !callback(&line) {
            // Kill the process if callback returns false
            let _ = child.kill();
            break;
        }
    }

    let status = child
        .wait()
        .map_err(|e| FastGAError::Other(format!("Failed to wait for FastGA: {e}")))?;

    if !status.success() {
        return Err(FastGAError::FastGAExecutionFailed(format!(
            "FastGA exited with status: {status}"
        )));
    }

    Ok(())
}
