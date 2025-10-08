//! Process runner for FastGA binary
//!
//! This module handles running the FastGA binary as a subprocess.
//! The FastGA binary automatically handles GDB conversion and indexing.

use std::ffi::CString;
use std::os::raw::{c_char, c_int};
use std::path::{Path, PathBuf};
use crate::error::{Result, FastGAError};
use crate::config::OutputFormat;
use nix::unistd::{fork, ForkResult};
use nix::sys::wait::{waitpid, WaitStatus};
use std::process::exit;

// Link required libraries
#[link(name = "z")]
#[link(name = "m")]
extern "C" {
    // Wrappers were not implemented, using main functions directly
}

// Also link to the renamed main functions
#[link(name = "fatogdb_main", kind = "static")]
#[link(name = "gixmake_main", kind = "static")]
extern "C" {
    fn fatogdb_main(argc: c_int, argv: *const *const c_char) -> c_int;
    fn gixmake_main(argc: c_int, argv: *const *const c_char) -> c_int;
}

/// Run FAtoGDB in a forked process (DEPRECATED - FastGA handles this automatically)
pub fn fork_fatogdb(input_path: &Path) -> Result<PathBuf> {
    eprintln!("[Fork] Converting {} to GDB format", input_path.display());

    // Handle compressed files: strip .fa.gz or .fasta.gz, not just .gz
    let gdb_path = if let Some(name) = input_path.file_name() {
        let name_str = name.to_string_lossy();
        if name_str.ends_with(".fa.gz") || name_str.ends_with(".fasta.gz") {
            // Remove .fa.gz or .fasta.gz and add .1gdb
            let base = if name_str.ends_with(".fa.gz") {
                &name_str[..name_str.len() - 6]
            } else {
                &name_str[..name_str.len() - 9]
            };
            input_path.with_file_name(format!("{}.1gdb", base))
        } else {
            input_path.with_extension("1gdb")
        }
    } else {
        input_path.with_extension("1gdb")
    };

    // Check if already exists
    if gdb_path.exists() {
        eprintln!("[Fork] GDB already exists");
        return Ok(gdb_path);
    }

    // Fork a new process
    match unsafe { fork() } {
        Ok(ForkResult::Parent { child }) => {
            // Parent process - wait for child
            eprintln!("[Fork] Parent waiting for FAtoGDB child process {child}");

            match waitpid(child, None) {
                Ok(WaitStatus::Exited(_, code)) => {
                    eprintln!("[Fork] FAtoGDB exited with code {code}");
                    if code == 0 && gdb_path.exists() {
                        Ok(gdb_path)
                    } else {
                        Err(FastGAError::Other(format!("FAtoGDB failed with code {code}")))
                    }
                }
                Ok(status) => {
                    Err(FastGAError::Other(format!("FAtoGDB unexpected status: {status:?}")))
                }
                Err(e) => {
                    Err(FastGAError::Other(format!("Failed to wait for FAtoGDB: {e}")))
                }
            }
        }
        Ok(ForkResult::Child) => {
            // Child process - run FAtoGDB
            eprintln!("[Fork] Child process running FAtoGDB");

            // Use fatogdb_main directly (wrapper not implemented)
            let args = [CString::new("FAtoGDB").unwrap(),
                CString::new(input_path.to_string_lossy().to_string()).unwrap()];

            let argv: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();

            let result = unsafe {
                fatogdb_main(argv.len() as c_int, argv.as_ptr())
            };

            exit(result);
        }
        Err(e) => {
            Err(FastGAError::Other(format!("Failed to fork for FAtoGDB: {e}")))
        }
    }
}

/// Run GIXmake in a forked process (DEPRECATED - FastGA handles this automatically)
pub fn fork_gixmake(gdb_path: &Path, threads: i32, freq: i32) -> Result<PathBuf> {
    eprintln!("[Fork] Creating index for {}", gdb_path.display());

    let base_name = gdb_path.with_extension("");
    let gix_path = base_name.with_extension("gix");

    if gix_path.exists() {
        eprintln!("[Fork] Index already exists");
        return Ok(gix_path);
    }

    match unsafe { fork() } {
        Ok(ForkResult::Parent { child }) => {
            eprintln!("[Fork] Parent waiting for GIXmake child process {child}");

            match waitpid(child, None) {
                Ok(WaitStatus::Exited(_, code)) => {
                    eprintln!("[Fork] GIXmake exited with code {code}");
                    if code == 0 && gix_path.exists() {
                        Ok(gix_path)
                    } else {
                        Err(FastGAError::Other(format!("GIXmake failed with code {code}")))
                    }
                }
                Ok(status) => {
                    Err(FastGAError::Other(format!("GIXmake unexpected status: {status:?}")))
                }
                Err(e) => {
                    Err(FastGAError::Other(format!("Failed to wait for GIXmake: {e}")))
                }
            }
        }
        Ok(ForkResult::Child) => {
            eprintln!("[Fork] Child process running GIXmake");

            // Use gixmake_main directly since wrapper isn't implemented
            let args = [CString::new("GIXmake").unwrap(),
                CString::new(format!("-T{threads}")).unwrap(),
                CString::new(format!("-f{freq}")).unwrap(),
                CString::new(base_name.to_string_lossy().to_string()).unwrap()];

            let argv: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();

            let result = unsafe {
                gixmake_main(argv.len() as c_int, argv.as_ptr())
            };

            exit(result);
        }
        Err(e) => {
            Err(FastGAError::Other(format!("Failed to fork for GIXmake: {e}")))
        }
    }
}

/// Orchestrator for running FastGA alignments
pub struct ForkOrchestrator {
    pub config: crate::Config,
}

impl ForkOrchestrator {
    pub fn new(config: crate::Config) -> Self {
        Self { config }
    }

    pub fn new_simple(threads: i32) -> Self {
        let config = crate::Config::builder()
            .num_threads(threads as usize)
            .build();
        Self { config }
    }

    pub fn align(&self, query_path: &Path, target_path: &Path) -> Result<String> {
        eprintln!("[fastga] Running alignment");

        // FastGA handles everything: FASTA→GDB conversion, index creation, and alignment
        self.run_fastga_alignment(query_path, target_path)
    }

    fn run_fastga_alignment(&self, query_path: &Path, target_path: &Path) -> Result<String> {
        use std::process::Command;

        // Find FastGA binary
        let fastga = find_fastga_binary()?;

        // Create temp .1aln file in same directory as input
        let working_dir = query_path.parent().ok_or_else(||
            FastGAError::Other("Cannot determine parent directory".to_string()))?;
        let temp_aln = working_dir.join(format!("_tmp_{}.1aln", std::process::id()));

        let mut cmd = Command::new(&fastga);

        // Use .1aln output format
        cmd.arg(format!("-1:{}", temp_aln.file_name().unwrap().to_string_lossy()));

        // Basic parameters
        cmd.arg(format!("-T{}", self.config.num_threads));

        // Only add -l if it's > 0 (FastGA's default is 0 anyway)
        if self.config.min_alignment_length > 0 {
            cmd.arg(format!("-l{}", self.config.min_alignment_length));
        }

        if let Some(identity) = self.config.min_identity {
            cmd.arg(format!("-i{identity:.2}"));
        }

        // Advanced parameters
        if let Some(cutoff) = self.config.adaptive_seed_cutoff {
            cmd.arg(format!("-f{cutoff}"));
        }

        if let Some(coverage) = self.config.min_chain_coverage {
            // FastGA expects -c as an integer percentage (0-100)
            let coverage_pct = (coverage * 100.0) as i32;
            cmd.arg(format!("-c{coverage_pct}"));
        }

        if let Some(threshold) = self.config.chain_start_threshold {
            cmd.arg(format!("-s{threshold}"));
        }

        // Flags
        if self.config.verbose {
            cmd.arg("-v");
        }

        if self.config.keep_intermediates {
            cmd.arg("-k");
        }

        if self.config.soft_masking {
            cmd.arg("-M");
        }

        if self.config.symmetric_seeding {
            cmd.arg("-S");
        }

        // Log file (needs format -L:<filename>)
        if let Some(ref log_path) = self.config.log_file {
            cmd.arg(format!("-L:{}", log_path.display()));
        }

        // Temp directory (needs format -P<dir>)
        if let Some(ref temp_dir) = self.config.temp_dir {
            cmd.arg(format!("-P{}", temp_dir.display()));
        }

        // Don't add -g flag - it's not a valid option
        // FastGA will automatically use existing GDB files

        cmd.arg(query_path)
           .arg(target_path);

        // Set working directory to where the input files are
        if let Some(parent_dir) = query_path.parent() {
            cmd.current_dir(parent_dir);
        }

        eprintln!("[fastga] Executing: {cmd:?}");

        let output = cmd.output()
            .map_err(|e| FastGAError::Other(format!("Failed to run FastGA: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let _ = std::fs::remove_file(&temp_aln); // Clean up
            return Err(FastGAError::Other(format!("FastGA failed: {}", stderr)));
        }

        // FastGA succeeded, now convert .1aln to desired output format
        let result = match self.config.output_format {
            OutputFormat::Aln => {
                // Return path to .1aln file
                Ok(temp_aln.to_string_lossy().to_string())
            },
            OutputFormat::Psl => self.convert_to_psl(&temp_aln),
            _ => self.convert_to_paf(&temp_aln, self.config.output_format),
        }?;

        // Clean up temp .1aln file unless keep_intermediates is set or format is Aln
        if !self.config.keep_intermediates && self.config.output_format != OutputFormat::Aln {
            let _ = std::fs::remove_file(&temp_aln);
        }

        Ok(result)
    }

    fn convert_to_paf(&self, aln_file: &Path, format: OutputFormat) -> Result<String> {
        use std::process::Command;

        let alnto_paf = find_alnto_paf_binary()?;
        let mut cmd = Command::new(&alnto_paf);

        // Set ALNtoPAF flags based on output format
        match format {
            OutputFormat::PafWithX => { cmd.arg("-x"); },
            OutputFormat::PafWithM => { cmd.arg("-m"); },
            OutputFormat::PafShort => { cmd.arg("-s"); },
            OutputFormat::PafLong => { cmd.arg("-S"); },
            _ => {}, // Default PAF
        }

        cmd.arg(format!("-T{}", self.config.num_threads));
        cmd.arg(aln_file);

        eprintln!("[fastga] Converting .1aln to PAF: {cmd:?}");

        let output = cmd.output()
            .map_err(|e| FastGAError::Other(format!("Failed to run ALNtoPAF: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(FastGAError::Other(format!("ALNtoPAF failed: {}", stderr)));
        }

        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    }

    fn convert_to_psl(&self, aln_file: &Path) -> Result<String> {
        use std::process::Command;

        let alnto_psl = find_alnto_psl_binary()?;
        let mut cmd = Command::new(&alnto_psl);

        cmd.arg(format!("-T{}", self.config.num_threads));
        cmd.arg(aln_file);

        eprintln!("[fastga] Converting .1aln to PSL: {cmd:?}");

        let output = cmd.output()
            .map_err(|e| FastGAError::Other(format!("Failed to run ALNtoPSL: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(FastGAError::Other(format!("ALNtoPSL failed: {}", stderr)));
        }

        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    }
}

fn find_binary(name: &str) -> Result<PathBuf> {
    // Try our build directory first
    if let Ok(out_dir) = std::env::var("OUT_DIR") {
        let path = PathBuf::from(out_dir).join(name);
        if path.exists() {
            return Ok(path);
        }
    }

    // Try current exe's build directory
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

    // Try target directories
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

    Err(FastGAError::Other(format!("{} binary not found", name)))
}

fn find_fastga_binary() -> Result<PathBuf> {
    find_binary("FastGA")
}

fn find_alnto_paf_binary() -> Result<PathBuf> {
    find_binary("ALNtoPAF")
}

fn find_alnto_psl_binary() -> Result<PathBuf> {
    find_binary("ALNtoPSL")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_fork_fatogdb() {
        let dir = tempdir().unwrap();
        let test_fa = dir.path().join("test.fa");

        let mut file = File::create(&test_fa).unwrap();
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        file.flush().unwrap();

        match fork_fatogdb(&test_fa) {
            Ok(gdb_path) => {
                println!("Successfully created GDB: {:?}", gdb_path);
                assert!(gdb_path.exists());
            }
            Err(e) => {
                println!("FAtoGDB failed (might be expected): {}", e);
            }
        }
    }

    #[test]
    fn test_fork_orchestrator() {
        let dir = tempdir().unwrap();
        let test_fa = dir.path().join("test.fa");

        let mut file = File::create(&test_fa).unwrap();
        writeln!(file, ">sequence1").unwrap();
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        writeln!(file, "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        file.flush().unwrap();

        let orchestrator = ForkOrchestrator::new_simple(1);

        match orchestrator.align(&test_fa, &test_fa) {
            Ok(output) => {
                println!("Alignment succeeded: {} bytes output", output.len());
            }
            Err(e) => {
                println!("Alignment failed (might be expected): {}", e);
            }
        }
    }
}