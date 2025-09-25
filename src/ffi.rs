//! FFI bindings to FastGA C library

use crate::error::{FastGAError, Result};
use std::ffi::CString;
use std::io::Read;
use std::os::raw::{c_char, c_int, c_void};
use std::os::unix::io::FromRawFd;
use std::path::Path;

// Link to the FastGA static libraries built by build.rs
#[link(name = "fastga_main", kind = "static")]
#[link(name = "fatogdb_main", kind = "static")]
#[link(name = "gixmake_main", kind = "static")]
#[link(name = "fastga_common", kind = "static")]
extern "C" {
    // FastGA's main entry point - we can call this directly
    fn fastga_main(argc: c_int, argv: *const *const c_char) -> c_int;

    // These would be the actual alignment functions if exposed
    // For now we'll use fastga_main as the entry point
}

/// Find the directory containing FastGA utility binaries
fn find_fastga_bin_dir() -> String {
    // First check environment variable
    if let Ok(dir) = std::env::var("FASTGA_BIN_DIR") {
        return dir;
    }

    // Try to find the OUT_DIR from the build
    // The binaries are in target/{debug,release}/build/fastga-rs-*/out/
    if let Ok(exe) = std::env::current_exe() {
        let mut path = exe.parent().unwrap().to_path_buf();

        // Go up to target dir
        while path.file_name() != Some(std::ffi::OsStr::new("target")) && path.parent().is_some() {
            path = path.parent().unwrap().to_path_buf();
        }

        // Find the build directory
        if path.file_name() == Some(std::ffi::OsStr::new("target")) {
            let profile = if cfg!(debug_assertions) {
                "debug"
            } else {
                "release"
            };
            let build_dir = path.join(profile).join("build");

            if build_dir.exists() {
                // Find fastga-rs-* directory
                if let Ok(entries) = std::fs::read_dir(&build_dir) {
                    for entry in entries {
                        if let Ok(entry) = entry {
                            let name = entry.file_name();
                            if name.to_string_lossy().starts_with("fastga-rs-") {
                                let out_dir = entry.path().join("out");
                                if out_dir.exists() && out_dir.join("FAtoGDB").exists() {
                                    return out_dir.to_string_lossy().to_string();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Fallback to current directory
    ".".to_string()
}

/// Safe wrapper for running FastGA alignment
pub fn run_fastga_alignment(
    query_path: &Path,
    target_path: &Path,
    min_identity: f64,
    min_length: i32,
    num_threads: i32,
) -> Result<Vec<u8>> {
    // Build argument array for FastGA
    let mut args = vec![
        CString::new("FastGA").unwrap(),
        CString::new("-pafx").unwrap(), // PAF with extended CIGAR
        CString::new(format!("-T{}", num_threads)).unwrap(),
        CString::new(format!("-l{}", min_length)).unwrap(),
    ];

    if min_identity > 0.0 {
        args.push(CString::new(format!("-i{:.2}", min_identity)).unwrap());
    }

    // Add file paths
    let query_str = query_path
        .to_str()
        .ok_or_else(|| FastGAError::Other("Invalid query path".to_string()))?;
    let target_str = target_path
        .to_str()
        .ok_or_else(|| FastGAError::Other("Invalid target path".to_string()))?;

    args.push(CString::new(query_str).unwrap());
    args.push(CString::new(target_str).unwrap());

    // Convert to raw pointers for C
    let c_args: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();

    // Set PATH to include our binary directory so FastGA can find FAtoGDB, GIXmake, etc.
    let bin_dir = find_fastga_bin_dir();

    // Save original PATH
    let original_path = std::env::var("PATH").unwrap_or_default();

    // Add our binary directory to PATH
    let new_path = format!("{}:{}", bin_dir, original_path);
    std::env::set_var("PATH", &new_path);

    // Create a pipe to capture stdout
    use std::io::Read;
    use std::os::unix::io::FromRawFd;

    let (read_pipe, write_pipe) = nix::unistd::pipe()
        .map_err(|e| FastGAError::Other(format!("Failed to create pipe: {}", e)))?;

    // Save original stdout
    let original_stdout = unsafe { libc::dup(1) };
    if original_stdout == -1 {
        return Err(FastGAError::Other("Failed to save stdout".to_string()));
    }

    // Redirect stdout to our pipe
    unsafe {
        libc::dup2(write_pipe, 1);
        libc::close(write_pipe);
    }

    // Call FastGA's main function
    let result = unsafe { fastga_main(c_args.len() as c_int, c_args.as_ptr()) };

    // Restore original stdout
    unsafe {
        libc::dup2(original_stdout, 1);
        libc::close(original_stdout);
    }

    // Restore original PATH
    std::env::set_var("PATH", original_path);

    if result != 0 {
        return Err(FastGAError::FastGAExecutionFailed(format!(
            "FastGA returned error code: {}",
            result
        )));
    }

    // Read output from pipe
    let mut output = Vec::new();
    let mut pipe_reader = unsafe { std::fs::File::from_raw_fd(read_pipe) };
    pipe_reader
        .read_to_end(&mut output)
        .map_err(|e| FastGAError::IoError(e))?;

    Ok(output)
}

/// Prepare a FASTA file as a GDB database
pub fn prepare_database(fasta_path: &Path) -> Result<String> {
    let gdb_path = fasta_path.with_extension("gdb");

    // Check if GDB already exists and is up to date
    if gdb_path.exists() {
        let fasta_modified = std::fs::metadata(fasta_path)?.modified()?;
        let gdb_modified = std::fs::metadata(&gdb_path)?.modified()?;

        if gdb_modified >= fasta_modified {
            // GDB is up to date
            return Ok(gdb_path.to_string_lossy().to_string());
        }
    }

    // Need to create GDB - for now we'll call FAtoGDB as a process
    // TODO: Call the C functions directly once we have proper bindings
    use std::process::Command;

    let output = Command::new(crate::embedded::get_binary_path("FAtoGDB")?)
        .arg(fasta_path)
        .output()?;

    if !output.status.success() {
        return Err(FastGAError::Other(format!(
            "Failed to create GDB: {}",
            String::from_utf8_lossy(&output.stderr)
        )));
    }

    Ok(gdb_path.to_string_lossy().to_string())
}
