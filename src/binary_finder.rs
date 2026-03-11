//! Shared utility for finding FastGA binaries
//!
//! This module provides a unified way to find FastGA utility binaries
//! that works in development, after cargo install, and with system installs.

use crate::error::{FastGAError, Result};
use std::path::PathBuf;

/// Walk up from a directory to find `target/`, then search
/// `target/{release,debug}/build/fastga-rs-*/out/` for the named binary.
fn find_in_target_build(start_dir: &std::path::Path, name: &str) -> Option<PathBuf> {
    let mut dir = start_dir.to_path_buf();
    // Walk up until we find a directory named "target"
    while dir.file_name().is_some_and(|n| n != "target") {
        if !dir.pop() {
            return None;
        }
    }
    if !dir.ends_with("target") {
        return None;
    }
    for profile in &["release", "debug"] {
        let build_dir = dir.join(profile).join("build");
        if let Ok(entries) = std::fs::read_dir(&build_dir) {
            for entry in entries.flatten() {
                if entry
                    .file_name()
                    .to_string_lossy()
                    .starts_with("fastga-rs-")
                {
                    let binary = entry.path().join("out").join(name);
                    if binary.exists() {
                        // Always return canonical (absolute) path
                        return Some(binary.canonicalize().unwrap_or(binary));
                    }
                }
            }
        }
    }
    None
}

/// Find a FastGA utility binary by name
///
/// Search order:
/// 1. Same directory as current executable (cargo install)
/// 2. OUT_DIR from build.rs (development only)
/// 3. Walk up from executable to target/ and scan build dirs
///
/// The system PATH is intentionally NOT searched to avoid version
/// mismatches with a globally installed FastGA.
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

    // 3. Walk up from current executable to find target/build/fastga-rs-*/out/
    if let Ok(exe_path) = std::env::current_exe() {
        let exe_path = exe_path.canonicalize().unwrap_or(exe_path);
        if let Some(exe_dir) = exe_path.parent() {
            if let Some(found) = find_in_target_build(exe_dir, name) {
                return Ok(found);
            }
        }
    }

    Err(FastGAError::Other(format!(
        "{name} binary not found. It should have been built by fastga-rs build.rs. \
         Try rebuilding with `cargo build --release`."
    )))
}

/// Get the directory containing FastGA binaries
///
/// This is used to set up PATH for FastGA's system() calls
pub fn get_binary_dir() -> Result<PathBuf> {
    let binary = find_binary("FastGA")?;
    binary
        .parent()
        .map(|p| p.to_path_buf())
        .ok_or_else(|| FastGAError::Other("Cannot determine FastGA binary directory".to_string()))
}
