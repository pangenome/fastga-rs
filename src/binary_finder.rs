//! Shared utility for finding FastGA binaries
//!
//! This module provides a unified way to find FastGA utility binaries
//! that works in development, after cargo install, and with system installs.

use crate::error::{FastGAError, Result};
use std::path::PathBuf;

/// Private install directory next to the running executable:
/// `<exe_dir>/../libexec/<exe_stem>/`.
///
/// Package managers use this to keep bundled FastGA binaries out of `bin/`,
/// where they would collide with a separately packaged FastGA.
fn libexec_dir(exe_path: &std::path::Path) -> Option<PathBuf> {
    let exe_dir = exe_path.parent()?;
    let stem = exe_path.file_stem()?;
    Some(exe_dir.parent()?.join("libexec").join(stem))
}

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
/// 2. `<exe_dir>/../libexec/<exe_stem>/` (packaged installs, e.g. conda)
/// 3. OUT_DIR from build.rs (development only)
/// 4. Walk up from executable to target/ and scan build dirs
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

        // 2. Private install directory, so packagers can keep these out of bin/
        if let Some(dir) = libexec_dir(&exe_path) {
            let binary = dir.join(name);
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

    // 4. Walk up from current executable to find target/build/fastga-rs-*/out/
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn libexec_dir_sits_beside_bin() {
        let got = libexec_dir(Path::new("/opt/conda/bin/impg")).unwrap();
        assert_eq!(got, PathBuf::from("/opt/conda/libexec/impg"));
    }

    #[test]
    fn libexec_dir_needs_a_parent_of_the_bin_dir() {
        assert_eq!(libexec_dir(Path::new("/impg")), None);
    }
}
