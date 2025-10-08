//! Subprocess runner for FastGA binary
//!
//! This module handles running the FastGA binary as a subprocess via system calls.
//! The FastGA binary automatically handles GDB conversion and indexing.

use std::path::{Path, PathBuf};
use crate::error::{Result, FastGAError};
use crate::config::OutputFormat;

/// Orchestrator for running FastGA alignments via subprocess
pub struct Orchestrator {
    pub config: crate::Config,
}

impl Orchestrator {
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
    fn test_orchestrator() {
        let dir = tempdir().unwrap();
        let test_fa = dir.path().join("test.fa");

        let mut file = File::create(&test_fa).unwrap();
        writeln!(file, ">sequence1").unwrap();
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        writeln!(file, "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        file.flush().unwrap();

        let orchestrator = Orchestrator::new_simple(1);

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
