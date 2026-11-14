//! Shared utility for finding FastGA binaries
//!
//! This module provides a unified way to find FastGA utility binaries
//! that works in development, after cargo install, and with system installs.

use crate::error::{FastGAError, Result};
use std::path::{Path, PathBuf};

/// Private install directory next to the running executable:
/// `<exe_dir>/../libexec/<exe_stem>/`.
///
/// Package managers use this to keep bundled FastGA binaries out of `bin/`,
/// where they would collide with a separately packaged FastGA.
fn libexec_dir(exe_path: &Path) -> Option<PathBuf> {
    let exe_dir = exe_path.parent()?;
    let stem = exe_path.file_stem()?;
    Some(exe_dir.parent()?.join("libexec").join(stem))
}

/// Newest `name` in the `out` directories of `krate` under `build`.
///
/// Covers both cargo layouts: `<krate>-<hash>/out` (cargo <= 1.99) and
/// `<krate>/<hash>/out` (cargo >= 1.100).
fn newest_in_out_dirs(build: &Path, krate: &str, name: &str) -> Option<PathBuf> {
    let subdirs = |dir: &Path| -> Vec<PathBuf> {
        std::fs::read_dir(dir)
            .map(|entries| entries.flatten().map(|e| e.path()).collect())
            .unwrap_or_default()
    };
    let old_prefix = format!("{krate}-");
    let mut hash_dirs = subdirs(&build.join(krate));
    hash_dirs.extend(subdirs(build).into_iter().filter(|p| {
        p.file_name()
            .is_some_and(|n| n.to_string_lossy().starts_with(&old_prefix))
    }));
    hash_dirs
        .into_iter()
        .map(|d| d.join("out").join(name))
        .filter_map(|bin| Some((bin.metadata().ok()?.modified().ok()?, bin)))
        .max_by_key(|(mtime, _)| *mtime)
        .map(|(_, bin)| bin)
}

/// Newest bundled `name` under one cargo `build` directory.
pub(crate) fn find_in_build_dir(build: &Path, name: &str) -> Option<PathBuf> {
    newest_in_out_dirs(build, "fastga-rs", name)
}

/// Walk up from a directory to find `target/`, then search
/// `target/{release,debug}/build/` for the named binary.
pub(crate) fn find_in_target_build(start_dir: &Path, name: &str) -> Option<PathBuf> {
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
    ["release", "debug"]
        .iter()
        .find_map(|profile| find_in_build_dir(&dir.join(profile).join("build"), name))
        // Always return canonical (absolute) path
        .map(|binary| binary.canonicalize().unwrap_or(binary))
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

    // 4. Walk up from current executable to target/, scan build dirs
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

    #[test]
    fn libexec_dir_sits_beside_bin() {
        let got = libexec_dir(Path::new("/opt/conda/bin/impg")).unwrap();
        assert_eq!(got, PathBuf::from("/opt/conda/libexec/impg"));
    }

    #[test]
    fn libexec_dir_needs_a_parent_of_the_bin_dir() {
        assert_eq!(libexec_dir(Path::new("/impg")), None);
    }

    #[test]
    fn finds_binary_in_both_cargo_layouts() {
        let tmp = tempfile::tempdir().unwrap();
        let build = tmp.path().join("build");
        let put = |rel: &str, secs: u64| {
            let bin = build.join(rel).join("out").join("FAtoGDB");
            std::fs::create_dir_all(bin.parent().unwrap()).unwrap();
            let file = std::fs::File::create(&bin).unwrap();
            let mtime = std::time::UNIX_EPOCH + std::time::Duration::from_secs(secs);
            file.set_modified(mtime).unwrap();
            bin
        };
        assert_eq!(find_in_build_dir(&build, "FAtoGDB"), None);

        // Decoys: other crates, no binary
        put("fastga-rs-sys/aaaa", 900);
        put("other-1111", 900);
        std::fs::create_dir_all(build.join("fastga-rs/eeee/out")).unwrap();
        assert_eq!(find_in_build_dir(&build, "FAtoGDB"), None);

        let old_layout = put("fastga-rs-1111", 100);
        assert_eq!(find_in_build_dir(&build, "FAtoGDB"), Some(old_layout));

        let new_layout = put("fastga-rs/2222", 200);
        assert_eq!(find_in_build_dir(&build, "FAtoGDB"), Some(new_layout));

        let newer_old = put("fastga-rs-3333", 300);
        assert_eq!(find_in_build_dir(&build, "FAtoGDB"), Some(newer_old));

        let newer_new = put("fastga-rs/4444", 400);
        assert_eq!(find_in_build_dir(&build, "FAtoGDB"), Some(newer_new));
        assert_eq!(find_in_build_dir(&build, "GIXmake"), None);
    }
}
