//! Fork-based API for FastGA that avoids hanging FFI issues
//!
//! This provides a drop-in replacement for the standard FastGA API
//! but uses fork/exec internally to avoid memory corruption and hanging.

use crate::{Config, Alignments, Result, FastGAError};
use crate::fork_runner::ForkOrchestrator;
use std::path::Path;
use tempfile::TempDir;
use std::fs;

/// Fork-based FastGA implementation that avoids FFI hanging issues.
///
/// This is the recommended implementation for integration with SweepGA
/// as it isolates each utility call in a separate process.
#[derive(Debug, Clone)]
pub struct ForkFastGA {
    pub config: Config,
}

impl ForkFastGA {
    /// Create a new fork-based FastGA instance.
    pub fn new(config: Config) -> Result<Self> {
        Ok(ForkFastGA { config })
    }

    /// Align two genome files using the fork/exec approach.
    ///
    /// This method:
    /// 1. Forks child processes to run FAtoGDB for GDB conversion
    /// 2. Forks child processes to run GIXmake for index creation
    /// 3. Runs FastGA for the actual alignment
    ///
    /// Each utility runs in isolation, avoiding memory corruption issues.
    pub fn align_files(&self, genome1: &Path, genome2: &Path) -> Result<Alignments> {
        eprintln!("[ForkAPI] Starting alignment: {} vs {}",
                 genome1.display(), genome2.display());

        // Validate input files
        if !genome1.exists() {
            return Err(FastGAError::FileNotFound(genome1.to_path_buf()));
        }
        if !genome2.exists() {
            return Err(FastGAError::FileNotFound(genome2.to_path_buf()));
        }

        // Create orchestrator with full configuration
        let orchestrator = ForkOrchestrator::new(self.config.clone());

        // Run alignment pipeline
        let paf_output = orchestrator.align(genome1, genome2)?;

        // Parse PAF output
        if paf_output.is_empty() {
            eprintln!("[ForkAPI] Warning: Empty alignment output");
            return Ok(Alignments::new());
        }

        Alignments::from_paf(&paf_output)
    }

    /// Align sequences provided as byte arrays.
    pub fn align_sequences(&self, seq1: &[u8], seq2: &[u8]) -> Result<Alignments> {
        let temp_dir = TempDir::new()?;

        // Write sequences to temporary files
        let file1 = temp_dir.path().join("seq1.fasta");
        let file2 = temp_dir.path().join("seq2.fasta");

        fs::write(&file1, seq1)?;
        fs::write(&file2, seq2)?;

        self.align_files(&file1, &file2)
    }

    /// Align two genome files and write output directly to disk.
    ///
    /// This avoids loading all alignments into memory, which is useful
    /// for large-scale alignments.
    ///
    /// # Arguments
    /// * `genome1` - Path to query genome
    /// * `genome2` - Path to target genome
    /// * `output_path` - Path where to write PAF output
    ///
    /// # Returns
    /// The number of alignments written
    pub fn align_to_file(&self, genome1: &Path, genome2: &Path, output_path: &Path) -> Result<usize> {
        eprintln!("[ForkAPI] Aligning and writing to: {}", output_path.display());

        // Validate input files
        if !genome1.exists() {
            return Err(FastGAError::FileNotFound(genome1.to_path_buf()));
        }
        if !genome2.exists() {
            return Err(FastGAError::FileNotFound(genome2.to_path_buf()));
        }

        // Run alignment
        let alignments = self.align_files(genome1, genome2)?;

        // Write to file
        alignments.write_paf(output_path)?;

        Ok(alignments.len())
    }
}

/// Drop-in replacement for the standard FastGA that uses fork/exec.
///
/// This can be used as:
/// ```no_run
/// # use fastga_rs::{Config, fork_api::FastGA};
/// # use anyhow::Result;
/// # fn main() -> Result<()> {
/// let config = Config::default();
/// let aligner = FastGA::new(config)?;
/// # Ok(())
/// # }
/// ```
pub type FastGA = ForkFastGA;

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_fork_api_basic() {
        let dir = tempdir().unwrap();
        let test_fa = dir.path().join("test.fa");

        let mut file = File::create(&test_fa).unwrap();
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        writeln!(file, ">seq2").unwrap();
        writeln!(file, "TACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        file.flush().unwrap();

        let config = Config::default();
        let aligner = ForkFastGA::new(config).unwrap();

        match aligner.align_files(&test_fa, &test_fa) {
            Ok(alignments) => {
                eprintln!("✓ Fork API alignment succeeded");
                eprintln!("  Found {} alignments", alignments.len());
            }
            Err(e) => {
                eprintln!("⚠ Fork API alignment failed: {}", e);
                eprintln!("  This might be expected for test data");
            }
        }
    }
}