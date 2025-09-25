//! Configuration options for FastGA alignment operations.
//!
//! This module provides a builder pattern for configuring alignment parameters,
//! allowing fine-tuned control over sensitivity, performance, and output.

use std::path::PathBuf;

/// Configuration for FastGA alignment operations.
///
/// This struct contains all parameters that control how FastGA performs alignments.
/// Use the builder pattern to construct configurations with non-default values.
///
/// # Default Values
/// - `min_alignment_length`: 100 bp
/// - `min_identity`: None (no filtering)
/// - `num_threads`: Number of CPU cores
/// - `chain_break`: 2000 (anti-diagonal distance)
/// - `chain_min`: 170
/// - `frequency`: 10 (k-mer frequency threshold)
#[derive(Debug, Clone)]
pub struct Config {
    /// Minimum alignment length in base pairs
    pub min_alignment_length: usize,

    /// Minimum identity fraction (0.0-1.0) for alignments
    pub min_identity: Option<f64>,

    /// Number of threads to use for alignment
    pub num_threads: usize,

    /// Maximum anti-diagonal distance for chaining
    pub chain_break: usize,

    /// Minimum chain score
    pub chain_min: usize,

    /// K-mer frequency threshold
    pub frequency: usize,

    /// Temporary directory for intermediate files
    pub temp_dir: Option<PathBuf>,

    /// Enable soft masking (lowercase sequences are masked)
    pub soft_masking: bool,

    /// Keep intermediate files for debugging
    pub keep_intermediates: bool,

    /// Enable verbose mode for detailed progress
    pub verbose: bool,

    /// Use symmetric seeding (default: false, not recommended)
    pub symmetric_seeding: bool,

    /// Log file path for detailed output
    pub log_file: Option<PathBuf>,

    /// Output format: "pafx" (default), "pafm", "pafs", "pafS", "psl"
    pub output_format: OutputFormat,

    /// Adaptive seed count cutoff (-f parameter)
    pub adaptive_seed_cutoff: Option<usize>,

    /// Minimum seed chain coverage in both genomes (-c parameter)
    pub min_chain_coverage: Option<f64>,

    /// Threshold for starting a new seed chain (-s parameter)
    pub chain_start_threshold: Option<usize>,
}

/// Output format for alignments
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OutputFormat {
    /// PAF with CIGAR string using X's for mismatches (default)
    PafWithX,
    /// PAF with CIGAR string using ='s for matches
    PafWithM,
    /// PAF with CS string in short form
    PafShort,
    /// PAF with CS string in long form
    PafLong,
    /// PSL format
    Psl,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            min_alignment_length: 100,
            min_identity: None,
            num_threads: num_cpus::get().max(1), // Ensure at least 1 thread
            chain_break: 2000,
            chain_min: 170,
            frequency: 10,
            temp_dir: None,
            soft_masking: true,
            keep_intermediates: false,
            verbose: false,
            symmetric_seeding: false,
            log_file: None,
            output_format: OutputFormat::PafWithX,
            adaptive_seed_cutoff: None,
            min_chain_coverage: None,
            chain_start_threshold: None,
        }
    }
}

impl Config {
    /// Creates a new configuration builder.
    ///
    /// # Example
    /// ```
    /// use fastga_rs::Config;
    ///
    /// let config = Config::builder()
    ///     .min_alignment_length(150)
    ///     .min_identity(0.8)
    ///     .num_threads(4)
    ///     .build();
    /// ```
    pub fn builder() -> ConfigBuilder {
        ConfigBuilder::default()
    }
}

/// Builder for constructing Config instances.
#[derive(Debug, Default)]
pub struct ConfigBuilder {
    config: Config,
}

impl ConfigBuilder {
    /// Sets the minimum alignment length.
    ///
    /// Alignments shorter than this will be filtered out.
    /// Default: 100 bp
    pub fn min_alignment_length(mut self, length: usize) -> Self {
        self.config.min_alignment_length = length;
        self
    }

    /// Sets the minimum identity fraction for alignments.
    ///
    /// Value should be between 0.0 and 1.0.
    /// Default: None (no identity filtering)
    pub fn min_identity(mut self, identity: f64) -> Self {
        assert!(
            identity >= 0.0 && identity <= 1.0,
            "Identity must be between 0.0 and 1.0"
        );
        self.config.min_identity = Some(identity);
        self
    }

    /// Sets the number of threads to use.
    ///
    /// Default: Number of CPU cores
    pub fn num_threads(mut self, threads: usize) -> Self {
        assert!(threads > 0, "Number of threads must be positive");
        self.config.num_threads = threads;
        self
    }

    /// Sets the chain break distance.
    ///
    /// This controls how far apart two alignments can be to still be chained together.
    /// Default: 2000
    pub fn chain_break(mut self, distance: usize) -> Self {
        self.config.chain_break = distance;
        self
    }

    /// Sets the minimum chain score.
    ///
    /// Chains with scores below this threshold are discarded.
    /// Default: 170
    pub fn chain_min(mut self, score: usize) -> Self {
        self.config.chain_min = score;
        self
    }

    /// Sets the k-mer frequency threshold.
    ///
    /// K-mers appearing more frequently than this are ignored.
    /// Default: 10
    pub fn frequency(mut self, freq: usize) -> Self {
        self.config.frequency = freq;
        self
    }

    /// Sets the temporary directory for intermediate files.
    ///
    /// Default: System temp directory
    pub fn temp_dir(mut self, path: PathBuf) -> Self {
        self.config.temp_dir = Some(path);
        self
    }

    /// Enables or disables soft masking.
    ///
    /// When enabled, lowercase sequences are treated as masked.
    /// Default: true
    pub fn soft_masking(mut self, enabled: bool) -> Self {
        self.config.soft_masking = enabled;
        self
    }

    /// Keep intermediate files for debugging.
    ///
    /// Default: false
    pub fn keep_intermediates(mut self, keep: bool) -> Self {
        self.config.keep_intermediates = keep;
        self
    }

    /// Enable verbose mode for detailed progress output.
    ///
    /// Default: false
    pub fn verbose(mut self, verbose: bool) -> Self {
        self.config.verbose = verbose;
        self
    }

    /// Use symmetric seeding (not recommended by FastGA authors).
    ///
    /// Default: false
    pub fn symmetric_seeding(mut self, symmetric: bool) -> Self {
        self.config.symmetric_seeding = symmetric;
        self
    }

    /// Set log file path for detailed output.
    ///
    /// Default: None
    pub fn log_file(mut self, path: PathBuf) -> Self {
        self.config.log_file = Some(path);
        self
    }

    /// Set output format for alignments.
    ///
    /// Default: PafWithX
    pub fn output_format(mut self, format: OutputFormat) -> Self {
        self.config.output_format = format;
        self
    }

    /// Set adaptive seed count cutoff (-f parameter).
    ///
    /// Default: None (use FastGA default)
    pub fn adaptive_seed_cutoff(mut self, cutoff: usize) -> Self {
        self.config.adaptive_seed_cutoff = Some(cutoff);
        self
    }

    /// Set minimum seed chain coverage in both genomes (-c parameter).
    ///
    /// Default: None (use FastGA default)
    pub fn min_chain_coverage(mut self, coverage: f64) -> Self {
        self.config.min_chain_coverage = Some(coverage);
        self
    }

    /// Set threshold for starting a new seed chain (-s parameter).
    ///
    /// Default: None (use FastGA default)
    pub fn chain_start_threshold(mut self, threshold: usize) -> Self {
        self.config.chain_start_threshold = Some(threshold);
        self
    }

    /// Builds the final Config instance.
    pub fn build(self) -> Config {
        self.config
    }
}

/// Preset configurations for common use cases.
impl Config {
    /// High-sensitivity configuration for distant homologs.
    ///
    /// - Lower minimum alignment length (50 bp)
    /// - No minimum identity threshold
    /// - Increased chain break distance
    pub fn high_sensitivity() -> Self {
        Config {
            min_alignment_length: 50,
            min_identity: None,
            chain_break: 3000,
            chain_min: 100,
            ..Default::default()
        }
    }

    /// Fast configuration for closely related genomes.
    ///
    /// - Higher minimum alignment length (200 bp)
    /// - Higher minimum identity (90%)
    /// - Stricter chaining parameters
    pub fn fast() -> Self {
        Config {
            min_alignment_length: 200,
            min_identity: Some(0.9),
            chain_break: 1000,
            chain_min: 250,
            frequency: 20,
            ..Default::default()
        }
    }

    /// Configuration optimized for repetitive genomes.
    ///
    /// - Higher frequency threshold to ignore common k-mers
    /// - Stricter chaining to avoid spurious alignments
    pub fn repetitive_genomes() -> Self {
        Config {
            frequency: 50,
            chain_break: 1500,
            chain_min: 200,
            ..Default::default()
        }
    }
}
