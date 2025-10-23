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
/// 2. $CARGO_HOME/lib/{parent_package}/ (for dependent packages after cargo install)
/// 3. OUT_DIR from build.rs (development only)
/// 4. Cargo build directories (development)
/// 5. System PATH
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

    // 2. Try $CARGO_HOME/lib/{parent_package}/ for packages that depend on fastga-rs
    // This allows dependent packages to install FastGA binaries to a known location
    let cargo_home = std::env::var("CARGO_HOME")
        .ok()
        .map(PathBuf::from)
        .or_else(|| {
            std::env::var("HOME")
                .ok()
                .map(|h| PathBuf::from(h).join(".cargo"))
        });

    if let Some(cargo_home) = cargo_home {
        let lib_dir = cargo_home.join("lib");
        // Check common package names that might use fastga-rs
        for package in &["sweepga", "fastga-rs", "fastga"] {
            let binary = lib_dir.join(package).join(name);
            if binary.exists() {
                return Ok(binary);
            }
        }
    }

    // 3. Try OUT_DIR (compile-time env var, only works during build)
    if let Ok(out_dir) = std::env::var("OUT_DIR") {
        let path = PathBuf::from(out_dir).join(name);
        if path.exists() {
            return Ok(path);
        }
    }

    // 4. Try current exe's build directory (development in target/)
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let build_dir = exe_dir.join("build");
            if let Ok(entries) = std::fs::read_dir(&build_dir) {
                for entry in entries.flatten() {
                    if entry
                        .file_name()
                        .to_string_lossy()
                        .starts_with("fastga-rs-")
                    {
                        let binary = entry.path().join("out").join(name);
                        if binary.exists() {
                            return Ok(binary);
                        }
                    }
                }
            }
        }
    }

    // 5. Try target directories (running from project root)
    for profile in &["debug", "release"] {
        let build_dir = PathBuf::from(format!("target/{profile}/build"));
        if let Ok(entries) = std::fs::read_dir(&build_dir) {
            for entry in entries.flatten() {
                if entry
                    .file_name()
                    .to_string_lossy()
                    .starts_with("fastga-rs-")
                {
                    let binary = entry.path().join("out").join(name);
                    if binary.exists() {
                        return Ok(binary);
                    }
                }
            }
        }
    }

    // 6. Fall back to PATH (system-installed FastGA)
    if let Ok(path) = which::which(name) {
        return Ok(path);
    }

    Err(FastGAError::Other(format!(
        "{name} binary not found. Install FastGA utilities or ensure they're in PATH."
    )))
}

/// Get the directory containing FastGA binaries
///
/// This is used to set up PATH for FastGA's system() calls
pub fn get_binary_dir() -> Result<PathBuf> {
    // Same search logic as find_binary, but return the directory

    // 1. Try same directory as the current executable
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let test_binary = exe_dir.join("FastGA");
            if test_binary.exists() {
                return Ok(exe_dir.to_path_buf());
            }
        }
    }

    // 2. Try $CARGO_HOME/lib/{parent_package}/
    let cargo_home = std::env::var("CARGO_HOME")
        .ok()
        .map(PathBuf::from)
        .or_else(|| {
            std::env::var("HOME")
                .ok()
                .map(|h| PathBuf::from(h).join(".cargo"))
        });

    if let Some(cargo_home) = cargo_home {
        let lib_dir = cargo_home.join("lib");
        for package in &["sweepga", "fastga-rs", "fastga"] {
            let package_lib = lib_dir.join(package);
            let test_binary = package_lib.join("FastGA");
            if test_binary.exists() {
                return Ok(package_lib);
            }
        }
    }

    // 3. Try OUT_DIR
    if let Ok(out_dir) = std::env::var("OUT_DIR") {
        let path = PathBuf::from(&out_dir);
        let test_binary = path.join("FastGA");
        if test_binary.exists() {
            return Ok(path);
        }
    }

    // 4. Try current exe's build directory
    if let Ok(exe_path) = std::env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let build_dir = exe_dir.join("build");
            if let Ok(entries) = std::fs::read_dir(&build_dir) {
                for entry in entries.flatten() {
                    if entry
                        .file_name()
                        .to_string_lossy()
                        .starts_with("fastga-rs-")
                    {
                        let out_dir = entry.path().join("out");
                        let test_binary = out_dir.join("FastGA");
                        if test_binary.exists() {
                            return Ok(out_dir);
                        }
                    }
                }
            }
        }
    }

    // 5. Try target directories
    for profile in &["debug", "release"] {
        let build_dir = PathBuf::from(format!("target/{profile}/build"));
        if let Ok(entries) = std::fs::read_dir(&build_dir) {
            for entry in entries.flatten() {
                if entry
                    .file_name()
                    .to_string_lossy()
                    .starts_with("fastga-rs-")
                {
                    let out_dir = entry.path().join("out");
                    let test_binary = out_dir.join("FastGA");
                    if test_binary.exists() {
                        return Ok(out_dir);
                    }
                }
            }
        }
    }

    Err(FastGAError::Other(
        "FastGA binary directory not found".to_string(),
    ))
}
