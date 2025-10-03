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
//! # use anyhow::Result;
//! # fn main() -> Result<()> {
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
//! # Ok(())
//! # }
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
pub mod aln_reader;
pub mod config;
pub mod embedded;
pub mod error;
pub mod ffi;
pub mod fork_api;
pub mod fork_runner;
pub mod intermediate;
pub mod onelib;
pub mod orchestrator;
pub mod query_set;
pub mod simple_runner;
pub mod streaming;
pub mod timeout;

use error::Result;
use std::path::Path;

pub use alignment::{Alignment, Alignments};
pub use aln_reader::{AlnReader, AlnRecord, AlnWriter};
pub use config::{Config, OutputFormat};
pub use error::FastGAError;
pub use query_set::{align_queries, QueryAlignmentIterator, QueryAlignmentSet};

/// Main interface to the FastGA alignment engine.
///
/// This uses the fork-based implementation to avoid FFI hanging issues.
/// Each utility (FAtoGDB, GIXmake) runs in an isolated child process.
#[derive(Debug)]
pub struct FastGA {
    inner: fork_api::ForkFastGA,
}

impl FastGA {
    /// Creates a new FastGA instance with the given configuration.
    ///
    /// # Arguments
    /// * `config` - Configuration parameters for alignment
    ///
    /// # Example
    /// ```no_run
    /// # use anyhow::Result;
    /// # fn main() -> Result<()> {
    /// use fastga_rs::{FastGA, Config};
    ///
    /// let config = Config::builder()
    ///     .min_alignment_length(150)
    ///     .min_identity(0.7)
    ///     .num_threads(8)
    ///     .build();
    ///
    /// let aligner = FastGA::new(config)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(config: Config) -> Result<Self> {
        Ok(FastGA {
            inner: fork_api::ForkFastGA::new(config)?,
        })
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
        self.inner.align_files(genome1, genome2)
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
        self.inner.align_sequences(seq1, seq2)
    }

    /// Align two genome files and write output directly to disk.
    ///
    /// # Arguments
    /// * `genome1` - Path to query genome
    /// * `genome2` - Path to target genome
    /// * `output_path` - Path where to write PAF output
    ///
    /// # Returns
    /// The number of alignments written
    pub fn align_to_file(&self, genome1: &Path, genome2: &Path, output_path: &Path) -> Result<usize> {
        self.inner.align_to_file(genome1, genome2, output_path)
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
