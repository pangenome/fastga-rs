//! FastGA binary management.
//!
//! This module handles the FastGA binaries that are built from source
//! during the cargo build process.

use crate::error::{FastGAError, Result};
use std::env;
use std::path::PathBuf;
use std::process::{Command, Stdio};

/// Get the path to a FastGA binary.
///
/// The binaries are built during cargo build and placed in OUT_DIR.
pub fn get_binary_path(binary_name: &str) -> Result<PathBuf> {
    // During build, binaries are in OUT_DIR
    if let Ok(out_dir) = env::var("OUT_DIR") {
        let path = PathBuf::from(out_dir).join(binary_name);
        if path.exists() {
            return Ok(path);
        }
    }

    // For development, also check target/debug and target/release
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());

    // Check debug build
    let debug_path = PathBuf::from(&manifest_dir)
        .join("target/debug/fastga_bins")
        .join(binary_name);
    if debug_path.exists() {
        return Ok(debug_path);
    }

    // Check release build
    let release_path = PathBuf::from(&manifest_dir)
        .join("target/release/fastga_bins")
        .join(binary_name);
    if release_path.exists() {
        return Ok(release_path);
    }

    // Last resort: check if it's in PATH
    if let Ok(output) = Command::new("which").arg(binary_name).output() {
        if output.status.success() {
            let path = String::from_utf8_lossy(&output.stdout).trim().to_string();
            return Ok(PathBuf::from(path));
        }
    }

    Err(FastGAError::Other(format!(
        "FastGA binary '{}' not found. Please run 'cargo build' first.",
        binary_name
    )))
}

/// Container for FastGA binary paths
pub struct FastGABinaries {
    pub fastga_path: PathBuf,
    pub alntopaf_path: PathBuf,
    pub fatogdb_path: PathBuf,
    pub gixmake_path: PathBuf,
}

impl FastGABinaries {
    /// Get paths to all FastGA binaries.
    pub fn new() -> Result<Self> {
        Ok(FastGABinaries {
            fastga_path: get_binary_path("FastGA")?,
            alntopaf_path: get_binary_path("ALNtoPAF")?,
            fatogdb_path: get_binary_path("FAtoGDB")?,
            gixmake_path: get_binary_path("GIXmake")?,
        })
    }

    /// Run FastGA with the given arguments.
    pub fn run_fastga(&self, args: &[&str]) -> Result<String> {
        let output = Command::new(&self.fastga_path)
            .args(args)
            .output()
            .map_err(|e| FastGAError::Other(format!("Failed to run FastGA: {}", e)))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(FastGAError::FastGAExecutionFailed(stderr.to_string()));
        }

        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    }

    /// Run FastGA with streaming output.
    pub fn run_fastga_streaming<F>(&self, args: &[&str], mut callback: F) -> Result<()>
    where
        F: FnMut(&str) -> bool,
    {
        use std::io::{BufRead, BufReader};

        let mut child = Command::new(&self.fastga_path)
            .args(args)
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .map_err(|e| FastGAError::Other(format!("Failed to spawn FastGA: {}", e)))?;

        let stdout = child
            .stdout
            .take()
            .ok_or_else(|| FastGAError::Other("Failed to capture stdout".to_string()))?;

        let reader = BufReader::new(stdout);

        for line in reader.lines() {
            let line = line.map_err(|e| FastGAError::IoError(e))?;
            if !callback(&line) {
                break;
            }
        }

        let status = child
            .wait()
            .map_err(|e| FastGAError::Other(format!("Failed to wait for FastGA: {}", e)))?;

        if !status.success() {
            return Err(FastGAError::FastGAExecutionFailed(format!(
                "FastGA exited with status: {}",
                status
            )));
        }

        Ok(())
    }
}

/// Run FastGA with the given arguments.
pub fn run_fastga(args: &[&str]) -> Result<String> {
    let fastga_path = get_binary_path("FastGA")?;

    let output = Command::new(&fastga_path)
        .args(args)
        .output()
        .map_err(|e| FastGAError::Other(format!("Failed to run FastGA: {}", e)))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FastGAError::FastGAExecutionFailed(stderr.to_string()));
    }

    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

/// Run FastGA with streaming output.
pub fn run_fastga_streaming<F>(args: &[&str], mut callback: F) -> Result<()>
where
    F: FnMut(&str) -> bool,
{
    use std::io::{BufRead, BufReader};

    let fastga_path = get_binary_path("FastGA")?;

    let mut child = Command::new(&fastga_path)
        .args(args)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| FastGAError::Other(format!("Failed to spawn FastGA: {}", e)))?;

    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| FastGAError::Other("Failed to capture stdout".to_string()))?;

    let reader = BufReader::new(stdout);

    for line in reader.lines() {
        let line = line.map_err(|e| FastGAError::IoError(e))?;
        if !callback(&line) {
            break;
        }
    }

    let status = child
        .wait()
        .map_err(|e| FastGAError::Other(format!("Failed to wait for FastGA: {}", e)))?;

    if !status.success() {
        return Err(FastGAError::FastGAExecutionFailed(format!(
            "FastGA exited with status: {}",
            status
        )));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test removed - embedded binaries approach no longer used
}
