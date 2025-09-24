//! Embedded FastGA binary management.
//!
//! This module handles the FastGA binaries that are compiled and embedded
//! directly into our Rust binary, providing a self-contained solution.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Once;
use tempfile::TempDir;
use crate::error::{Result, FastGAError};

// Include the compiled FastGA binaries as bytes
const FASTGA_BINARY: &[u8] = include_bytes!(env!("FASTGA_BINARY"));
const ALNTOPAF_BINARY: &[u8] = include_bytes!(env!("ALNTOPAF_BINARY"));
const FATOGDB_BINARY: &[u8] = include_bytes!(env!("FATOGDB_BINARY"));

static INIT: Once = Once::new();
static mut BINARY_DIR: Option<PathBuf> = None;

/// Container for extracted FastGA binaries
pub struct EmbeddedFastGA {
    temp_dir: Option<TempDir>,
    fastga_path: PathBuf,
    alntopaf_path: PathBuf,
    fatogdb_path: PathBuf,
}

impl EmbeddedFastGA {
    /// Initialize embedded FastGA binaries.
    ///
    /// This extracts the embedded binaries to a temporary directory
    /// and makes them executable.
    pub fn new() -> Result<Self> {
        // Create temporary directory for binaries
        let temp_dir = TempDir::new()
            .map_err(|e| FastGAError::Other(format!("Failed to create temp dir: {}", e)))?;

        let dir_path = temp_dir.path();

        // Extract FastGA binary
        let fastga_path = dir_path.join("FastGA");
        Self::extract_binary(FASTGA_BINARY, &fastga_path)?;

        // Extract ALNtoPAF binary
        let alntopaf_path = dir_path.join("ALNtoPAF");
        Self::extract_binary(ALNTOPAF_BINARY, &alntopaf_path)?;

        // Extract FAtoGDB binary
        let fatogdb_path = dir_path.join("FAtoGDB");
        Self::extract_binary(FATOGDB_BINARY, &fatogdb_path)?;

        Ok(EmbeddedFastGA {
            temp_dir: Some(temp_dir),
            fastga_path,
            alntopaf_path,
            fatogdb_path,
        })
    }

    /// Extract a binary from embedded bytes to a file.
    fn extract_binary(binary_data: &[u8], path: &Path) -> Result<()> {
        let mut file = fs::File::create(path)
            .map_err(|e| FastGAError::Other(format!("Failed to create binary file: {}", e)))?;

        file.write_all(binary_data)
            .map_err(|e| FastGAError::Other(format!("Failed to write binary: {}", e)))?;

        // Make executable on Unix
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perms = fs::metadata(path)?.permissions();
            perms.set_mode(0o755);
            fs::set_permissions(path, perms)?;
        }

        Ok(())
    }

    /// Get path to the FastGA binary.
    pub fn fastga_path(&self) -> &Path {
        &self.fastga_path
    }

    /// Get path to the ALNtoPAF binary.
    pub fn alntopaf_path(&self) -> &Path {
        &self.alntopaf_path
    }

    /// Get path to the FAtoGDB binary.
    pub fn fatogdb_path(&self) -> &Path {
        &self.fatogdb_path
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

        let stdout = child.stdout.take()
            .ok_or_else(|| FastGAError::Other("Failed to capture stdout".to_string()))?;

        let reader = BufReader::new(stdout);

        for line in reader.lines() {
            let line = line.map_err(|e| FastGAError::IoError(e))?;
            if !callback(&line) {
                break;
            }
        }

        let status = child.wait()
            .map_err(|e| FastGAError::Other(format!("Failed to wait for FastGA: {}", e)))?;

        if !status.success() {
            return Err(FastGAError::FastGAExecutionFailed(
                format!("FastGA exited with status: {}", status)
            ));
        }

        Ok(())
    }
}

/// Global singleton instance of embedded FastGA.
///
/// This ensures we only extract the binaries once per process.
pub fn get_embedded_fastga() -> Result<&'static EmbeddedFastGA> {
    static mut INSTANCE: Option<EmbeddedFastGA> = None;

    unsafe {
        INIT.call_once(|| {
            match EmbeddedFastGA::new() {
                Ok(fastga) => {
                    // Leak the temp directory so it lives for the entire program
                    let mut leaked = fastga;
                    let temp_dir = leaked.temp_dir.take();
                    if let Some(dir) = temp_dir {
                        // This prevents the temp directory from being deleted
                        std::mem::forget(dir);
                    }
                    INSTANCE = Some(leaked);
                }
                Err(e) => {
                    eprintln!("Failed to initialize embedded FastGA: {}", e);
                }
            }
        });

        INSTANCE.as_ref()
            .ok_or_else(|| FastGAError::Other("Failed to initialize embedded FastGA".to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_embedded_extraction() {
        let fastga = EmbeddedFastGA::new().unwrap();
        assert!(fastga.fastga_path().exists());
        assert!(fastga.alntopaf_path().exists());
        assert!(fastga.fatogdb_path().exists());
    }

    #[test]
    fn test_embedded_singleton() {
        let fastga1 = get_embedded_fastga().unwrap();
        let fastga2 = get_embedded_fastga().unwrap();
        // Should be the same instance
        assert_eq!(fastga1 as *const _, fastga2 as *const _);
    }
}