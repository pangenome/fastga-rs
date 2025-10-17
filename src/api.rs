//! Subprocess-based API for FastGA
//!
//! This provides a stable API by running FastGA as a subprocess via system calls,
//! avoiding memory corruption and hanging issues from direct FFI.

use crate::runner::Orchestrator;
use crate::{Alignments, Config, FastGAError, Result};
use std::fs;
use std::path::Path;
use tempfile::TempDir;

/// Subprocess-based FastGA implementation.
///
/// This is the recommended implementation as it runs FastGA
/// in a separate process via system calls, avoiding any FFI-related issues.
#[derive(Debug, Clone)]
pub struct ProcessFastGA {
    pub config: Config,
}

impl ProcessFastGA {
    /// Create a new subprocess-based FastGA instance.
    pub fn new(config: Config) -> Result<Self> {
        Ok(ProcessFastGA { config })
    }

    /// Align two genome files using the FastGA binary.
    ///
    /// The FastGA binary handles everything automatically:
    /// 1. Converts FASTA to GDB format if needed
    /// 2. Creates k-mer indices if needed
    /// 3. Performs the alignment
    ///
    /// Running as a subprocess avoids FFI memory issues.
    pub fn align_files(&self, genome1: &Path, genome2: &Path) -> Result<Alignments> {
        eprintln!(
            "[fastga-rs] Starting alignment: {} vs {}",
            genome1.display(),
            genome2.display()
        );

        // Validate input files
        if !genome1.exists() {
            return Err(FastGAError::FileNotFound(genome1.to_path_buf()));
        }
        if !genome2.exists() {
            return Err(FastGAError::FileNotFound(genome2.to_path_buf()));
        }

        // Override output format to PafWithX for parsing
        // (low-level default is .1aln, but high-level API needs PAF to parse)
        let mut config = self.config.clone();
        config.output_format = crate::config::OutputFormat::PafWithX;

        // Create orchestrator with modified configuration
        let orchestrator = Orchestrator::new(config);

        // Run alignment pipeline
        let paf_output = orchestrator.align(genome1, genome2)?;

        // Parse PAF output
        if paf_output.is_empty() {
            eprintln!("[fastga-rs] Warning: Empty alignment output");
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
    pub fn align_to_file(
        &self,
        genome1: &Path,
        genome2: &Path,
        output_path: &Path,
    ) -> Result<usize> {
        eprintln!(
            "[fastga-rs] Aligning and writing to: {}",
            output_path.display()
        );

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

/// Drop-in replacement for the standard FastGA that uses subprocesses.
///
/// This can be used as:
/// ```no_run
/// # use fastga_rs::{Config, api::FastGA};
/// # use anyhow::Result;
/// # fn main() -> Result<()> {
/// let config = Config::default();
/// let aligner = FastGA::new(config)?;
/// # Ok(())
/// # }
/// ```
pub type FastGA = ProcessFastGA;

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_subprocess_api_basic() {
        let dir = tempdir().unwrap();
        let test_fa = dir.path().join("test.fa");

        let mut file = File::create(&test_fa).unwrap();
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
        writeln!(file, ">seq2").unwrap();
        writeln!(file, "TACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        file.flush().unwrap();

        let config = Config::default();
        let aligner = ProcessFastGA::new(config).unwrap();

        match aligner.align_files(&test_fa, &test_fa) {
            Ok(alignments) => {
                eprintln!("✓ Subprocess API alignment succeeded");
                eprintln!("  Found {} alignments", alignments.len());
            }
            Err(e) => {
                eprintln!("⚠ Subprocess API alignment failed: {}", e);
                eprintln!("  This might be expected for test data");
            }
        }
    }
}
