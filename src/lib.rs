//! # FastGA-RS: Rust Bindings for FastGA Genome Aligner
//!
//! This library provides safe Rust bindings for the FastGA genome alignment tool,
//! enabling efficient local alignment of high-quality genome assemblies.
//!
//! ## Overview
//!
//! FastGA-RS allows you to:
//! - Perform fast local DNA alignments between two genomes
//! - Generate alignments in PAF format with extended CIGAR strings (using '=' and 'X')
//! - Handle large genome assemblies efficiently
//! - Configure alignment parameters for different sensitivity levels
//!
//! ## Key Features
//!
//! - **Extended CIGAR Support**: Produces detailed alignment strings with explicit
//!   match ('=') and mismatch ('X') operators, crucial for accurate variant calling
//! - **High Performance**: Leverages FastGA's novel adaptive seed finding algorithm
//! - **Memory Efficient**: Uses trace point compression for alignment storage
//! - **Thread-Safe**: Supports parallel alignment operations
//!
//! ## Example Usage
//!
//! ```no_run
//! use fastga_rs::{FastGA, Config};
//! use std::path::Path;
//!
//! // Create aligner with default configuration
//! let aligner = FastGA::new(Config::default())?;
//!
//! // Align two genome files
//! let alignments = aligner.align_files(
//!     Path::new("genome1.fasta"),
//!     Path::new("genome2.fasta")
//! )?;
//!
//! // Get alignments in PAF format with extended CIGAR
//! let paf_output = alignments.to_paf()?;
//! println!("{}", paf_output);
//! ```
//!
//! ## Architecture
//!
//! The library is structured in several modules:
//! - `ffi`: Low-level FFI bindings to FastGA C code
//! - `config`: Configuration options for alignment parameters
//! - `alignment`: High-level alignment representation and format conversion
//! - `error`: Error types for the library
//!
//! ## Thread Safety
//!
//! The FastGA aligner is thread-safe for read operations but alignment operations
//! should be synchronized when using the same aligner instance from multiple threads.

pub mod alignment;
pub mod config;
pub mod error;
mod ffi;

use std::path::Path;
use std::fs;
use std::process::Command;
use tempfile::TempDir;
use error::Result;

pub use alignment::{Alignment, Alignments};
pub use config::Config;
pub use error::FastGAError;

/// Main interface to the FastGA alignment engine.
///
/// This struct provides the primary API for performing genome alignments.
/// It manages temporary files, executes the FastGA binary with appropriate
/// parameters, and parses the output into structured alignment data.
#[derive(Debug)]
pub struct FastGA {
    config: Config,
}

impl FastGA {
    /// Creates a new FastGA instance with the given configuration.
    ///
    /// # Arguments
    /// * `config` - Configuration parameters for alignment
    ///
    /// # Example
    /// ```no_run
    /// use fastga_rs::{FastGA, Config};
    ///
    /// let config = Config::builder()
    ///     .min_alignment_length(150)
    ///     .min_identity(0.7)
    ///     .num_threads(8)
    ///     .build();
    ///
    /// let aligner = FastGA::new(config)?;
    /// ```
    pub fn new(config: Config) -> Result<Self> {
        Ok(FastGA { config })
    }

    /// Aligns two genome files and returns the alignments.
    ///
    /// This method takes paths to two FASTA files, runs FastGA alignment,
    /// and returns structured alignment data with extended CIGAR strings.
    ///
    /// # Arguments
    /// * `genome1` - Path to the first genome file (FASTA format)
    /// * `genome2` - Path to the second genome file (FASTA format)
    ///
    /// # Returns
    /// A collection of alignments between the two genomes
    ///
    /// # Errors
    /// Returns an error if:
    /// - Input files don't exist or aren't readable
    /// - FastGA execution fails
    /// - Output parsing fails
    pub fn align_files(&self, genome1: &Path, genome2: &Path) -> Result<Alignments> {
        // Validate input files
        if !genome1.exists() {
            return Err(FastGAError::FileNotFound(genome1.to_path_buf()));
        }
        if !genome2.exists() {
            return Err(FastGAError::FileNotFound(genome2.to_path_buf()));
        }

        // Create temporary directory for intermediate files
        let temp_dir = TempDir::new()?;

        // Execute FastGA with extended CIGAR output
        let output = self.run_fastga(genome1, genome2, &temp_dir)?;

        // Parse PAF output into alignments
        Alignments::from_paf(&output)
    }

    /// Aligns sequences provided as byte arrays.
    ///
    /// This method creates temporary FASTA files from the provided sequences,
    /// performs alignment, and returns the results.
    ///
    /// # Arguments
    /// * `seq1` - First sequence as bytes
    /// * `seq2` - Second sequence as bytes
    ///
    /// # Returns
    /// Alignments between the two sequences
    pub fn align_sequences(&self, seq1: &[u8], seq2: &[u8]) -> Result<Alignments> {
        let temp_dir = TempDir::new()?;

        // Write sequences to temporary files
        let file1 = temp_dir.path().join("seq1.fasta");
        let file2 = temp_dir.path().join("seq2.fasta");

        fs::write(&file1, seq1)?;
        fs::write(&file2, seq2)?;

        self.align_files(&file1, &file2)
    }

    /// Executes FastGA binary with appropriate parameters.
    fn run_fastga(&self, genome1: &Path, genome2: &Path, _temp_dir: &TempDir) -> Result<String> {
        // Build FastGA command with -pafx flags for PAF output with extended CIGAR
        let mut cmd = Command::new("FastGA");

        cmd.arg("-pafx")  // PAF output with extended CIGAR
            .arg("-T").arg(self.config.num_threads.to_string())
            .arg("-l").arg(self.config.min_alignment_length.to_string());

        if let Some(min_id) = self.config.min_identity {
            cmd.arg("-e").arg(format!("{:.2}", 1.0 - min_id));
        }

        cmd.arg(genome1)
            .arg(genome2);

        let output = cmd.output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(FastGAError::FastGAExecutionFailed(stderr.to_string()));
        }

        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    }

    /// Aligns sequences with streaming output and filtering capability.
    ///
    /// This method provides fine-grained control over alignment processing,
    /// allowing you to filter or process alignments as they are generated.
    /// This is particularly useful for:
    /// - Large-scale alignments where memory is a concern
    /// - Query-vs-all scenarios where you want to filter per query
    /// - Real-time alignment processing pipelines
    ///
    /// # Arguments
    /// * `genome1` - Path to query genome file
    /// * `genome2` - Path to target genome file
    /// * `callback` - Function called for each alignment as it's generated
    ///
    /// # Example
    /// ```no_run
    /// use fastga_rs::{FastGA, Config, Alignment};
    ///
    /// let aligner = FastGA::new(Config::default())?;
    /// let mut filtered_alignments = Vec::new();
    ///
    /// aligner.align_streaming(
    ///     Path::new("query.fasta"),
    ///     Path::new("target.fasta"),
    ///     |alignment: Alignment| {
    ///         // Filter alignments with identity > 90%
    ///         if alignment.identity() > 0.9 {
    ///             filtered_alignments.push(alignment);
    ///             true  // Continue processing
    ///         } else {
    ///             true  // Continue but don't store
    ///         }
    ///     }
    /// )?;
    /// ```
    pub fn align_streaming<F>(&self, genome1: &Path, genome2: &Path, mut callback: F) -> Result<()>
    where
        F: FnMut(Alignment) -> bool,
    {
        // This will use the FFI bindings to hook into FastGA's internal
        // alignment generation and call the callback for each alignment
        // TODO: Implement using FFI bindings
        unimplemented!("Streaming alignment API will be implemented with FFI bindings")
    }

    /// Aligns a single query sequence against all sequences in a target file.
    ///
    /// This method is optimized for the common case of aligning one query
    /// against a database of targets, with immediate filtering capability.
    ///
    /// # Arguments
    /// * `query` - Query sequence as bytes
    /// * `target_file` - Path to target genome file
    /// * `filter` - Optional filter function for alignments
    ///
    /// # Returns
    /// Only the alignments that pass the filter
    pub fn align_query_vs_all<F>(&self, query: &[u8], target_file: &Path, filter: Option<F>) -> Result<Vec<Alignment>>
    where
        F: Fn(&Alignment) -> bool,
    {
        // TODO: Implement using lower-level FFI bindings
        unimplemented!("Query-vs-all API will be implemented with FFI bindings")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_builder() {
        let config = Config::builder()
            .min_alignment_length(200)
            .min_identity(0.8)
            .num_threads(4)
            .build();

        assert_eq!(config.min_alignment_length, 200);
        assert_eq!(config.min_identity, Some(0.8));
        assert_eq!(config.num_threads, 4);
    }
}