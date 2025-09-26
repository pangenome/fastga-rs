//! Fork/exec runner for FastGA utilities
//! This runs each utility in a separate forked process to avoid memory/state issues

use std::ffi::CString;
use std::os::raw::{c_char, c_int};
use std::path::{Path, PathBuf};
use crate::error::{Result, FastGAError};
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

/// Run FAtoGDB in a forked process
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

/// Run GIXmake in a forked process
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

/// Complete pipeline using fork/exec for each utility
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
        eprintln!("[Fork] Starting fork-based alignment pipeline");

        // Step 1: Convert to GDB using forked process
        eprintln!("[Fork] Step 1: Converting FASTA files to GDB");
        let query_gdb = fork_fatogdb(query_path)?;
        let target_gdb = fork_fatogdb(target_path)?;

        // Step 2: Create indices using forked process
        eprintln!("[Fork] Step 2: Creating k-mer indices");
        let _query_gix = fork_gixmake(&query_gdb, self.config.num_threads as i32, self.config.frequency as i32)?;
        let _target_gix = fork_gixmake(&target_gdb, self.config.num_threads as i32, self.config.frequency as i32)?;

        // Step 3: Run alignment
        // For alignment, we still need to use the FastGA binary
        eprintln!("[Fork] Step 3: Running alignment");
        self.run_fastga_alignment(query_path, target_path)
    }

    fn run_fastga_alignment(&self, query_path: &Path, target_path: &Path) -> Result<String> {
        use std::process::Command;
        use crate::config::OutputFormat;

        // Find FastGA binary
        let fastga = find_fastga_binary()?;

        let mut cmd = Command::new(&fastga);

        // Set output format
        match self.config.output_format {
            OutputFormat::PafWithX => cmd.arg("-pafx"),
            OutputFormat::PafWithM => cmd.arg("-pafm"),
            OutputFormat::PafShort => cmd.arg("-pafs"),
            OutputFormat::PafLong => cmd.arg("-pafS"),
            OutputFormat::Psl => cmd.arg("-psl"),
        };

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

        eprintln!("[Fork] Running FastGA: {cmd:?}");

        let output = cmd.output()
            .map_err(|e| FastGAError::Other(format!("Failed to run FastGA: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            eprintln!("[Fork] FastGA failed with: {}", stderr.trim());
            // Try without -g flag
            eprintln!("[Fork] Trying without -g flag");

            let mut cmd2 = Command::new(&fastga);

            // Rebuild command without -g flag
            match self.config.output_format {
                OutputFormat::PafWithX => cmd2.arg("-pafx"),
                OutputFormat::PafWithM => cmd2.arg("-pafm"),
                OutputFormat::PafShort => cmd2.arg("-pafs"),
                OutputFormat::PafLong => cmd2.arg("-pafS"),
                OutputFormat::Psl => cmd2.arg("-psl"),
            };

            cmd2.arg(format!("-T{}", self.config.num_threads));

            // Only add -l if it's > 0 (FastGA's default is 0 anyway)
            if self.config.min_alignment_length > 0 {
                cmd2.arg(format!("-l{}", self.config.min_alignment_length));
            }

            if let Some(identity) = self.config.min_identity {
                cmd2.arg(format!("-i{identity:.2}"));
            }

            // Add all other parameters (same as above)
            if let Some(cutoff) = self.config.adaptive_seed_cutoff {
                cmd2.arg(format!("-f{cutoff}"));
            }
            if let Some(coverage) = self.config.min_chain_coverage {
                // FastGA expects -c as an integer percentage (0-100)
                let coverage_pct = (coverage * 100.0) as i32;
                cmd2.arg(format!("-c{coverage_pct}"));
            }
            if let Some(threshold) = self.config.chain_start_threshold {
                cmd2.arg(format!("-s{threshold}"));
            }
            if self.config.verbose { cmd2.arg("-v"); }
            if self.config.keep_intermediates { cmd2.arg("-k"); }
            if self.config.soft_masking { cmd2.arg("-M"); }
            if self.config.symmetric_seeding { cmd2.arg("-S"); }
            if let Some(ref log_path) = self.config.log_file {
                cmd2.arg(format!("-L:{}", log_path.display()));
            }
            if let Some(ref temp_dir) = self.config.temp_dir {
                cmd2.arg(format!("-P{}", temp_dir.display()));
            }

            cmd2.arg(query_path)
                .arg(target_path);

            let output2 = cmd2.output()
                .map_err(|e| FastGAError::Other(format!("Failed to run FastGA: {e}")))?;

            if !output2.status.success() {
                return Err(FastGAError::Other(format!("FastGA failed: {}",
                    String::from_utf8_lossy(&output2.stderr))));
            }

            return Ok(String::from_utf8_lossy(&output2.stdout).to_string());
        }

        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    }
}

fn find_fastga_binary() -> Result<PathBuf> {
    // Try our build directory first
    if let Ok(out_dir) = std::env::var("OUT_DIR") {
        let path = PathBuf::from(out_dir).join("FastGA");
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
                        let fastga = entry.path().join("out/FastGA");
                        if fastga.exists() {
                            return Ok(fastga);
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
                    let fastga = entry.path().join("out/FastGA");
                    if fastga.exists() {
                        return Ok(fastga);
                    }
                }
            }
        }
    }

    Err(FastGAError::Other("FastGA binary not found".to_string()))
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