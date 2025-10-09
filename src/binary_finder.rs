//! Shared utility for finding FastGA binaries
//!
//! This module provides a unified way to find FastGA utility binaries
//! that works in development, after cargo install, and with system installs.

use crate::error::{FastGAError, Result};
use std::path::PathBuf;

/// Find a FastGA utility binary by name
///
/// Search order:
/// 1. Same directory as current executable (cargo install)
/// 2. OUT_DIR from build.rs (development only)
/// 3. Cargo build directories (development)
/// 4. System PATH
pub fn find_binary(name: &str) -> Result<PathBuf> {
    // 1. Try same directory as the current executable (for cargo install)
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let binary = exe_dir.join(name);
            if binary.exists() {
                return Ok(binary);
            }
        }
    }

    // 2. Try OUT_DIR (compile-time env var, only works during build)
    if let Ok(out_dir) = std::env::var("OUT_DIR") {
        let path = PathBuf::from(out_dir).join(name);
        if path.exists() {
            return Ok(path);
        }
    }

    // 3. Try current exe's build directory (development in target/)
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let build_dir = exe_dir.join("build");
            if let Ok(entries) = std::fs::read_dir(&build_dir) {
                for entry in entries.flatten() {
                    if entry.file_name().to_string_lossy().starts_with("fastga-rs-") {
                        let binary = entry.path().join("out").join(name);
                        if binary.exists() {
                            return Ok(binary);
                        }
                    }
                }
            }
        }
    }

    // 4. Try target directories (running from project root)
    for profile in &["debug", "release"] {
        let build_dir = PathBuf::from(format!("target/{profile}/build"));
        if let Ok(entries) = std::fs::read_dir(&build_dir) {
            for entry in entries.flatten() {
                if entry.file_name().to_string_lossy().starts_with("fastga-rs-") {
                    let binary = entry.path().join("out").join(name);
                    if binary.exists() {
                        return Ok(binary);
                    }
                }
            }
        }
    }

    // 5. Fall back to PATH (system-installed FastGA)
    if let Ok(path) = which::which(name) {
        return Ok(path);
    }

    Err(FastGAError::Other(format!(
        "{} binary not found. Install FastGA utilities or ensure they're in PATH.",
        name
    )))
}
