//! API for exposing intermediate steps in the alignment pipeline
//! This allows for debugging and fine-grained control over the alignment process

use crate::error::{FastGAError, Result};
use crate::Config;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Progress callback type for intermediate pipeline steps
type ProgressCallback = Box<dyn Fn(&str, &str)>;

/// Intermediate pipeline steps for alignment
pub struct AlignmentPipeline {
    config: Config,
    progress: Option<ProgressCallback>,
}

impl AlignmentPipeline {
    pub fn new(config: Config) -> Self {
        Self {
            config,
            progress: None,
        }
    }

    /// Add progress callback
    pub fn with_progress<F>(mut self, callback: F) -> Self
    where
        F: Fn(&str, &str) + 'static,
    {
        self.progress = Some(Box::new(callback));
        self
    }

    fn report_progress(&self, stage: &str, message: &str) {
        if let Some(ref callback) = self.progress {
            callback(stage, message);
        }
        eprintln!("[Pipeline] {stage}: {message}");
    }

    /// Step 1: Validate input files
    pub fn validate_inputs(&self, query: &Path, target: &Path) -> Result<()> {
        self.report_progress("validate", "Checking input files");

        if !query.exists() {
            return Err(FastGAError::FileNotFound(query.to_path_buf()));
        }
        if !target.exists() {
            return Err(FastGAError::FileNotFound(target.to_path_buf()));
        }

        // Check file size
        let query_size = std::fs::metadata(query)?.len();
        let target_size = std::fs::metadata(target)?.len();

        self.report_progress(
            "validate",
            &format!("Query: {query_size} bytes, Target: {target_size} bytes"),
        );

        if query_size == 0 || target_size == 0 {
            return Err(FastGAError::Other(
                "Input files cannot be empty".to_string(),
            ));
        }

        Ok(())
    }

    /// Step 2: Prepare database files (.gdb format)
    pub fn prepare_database(&self, fasta_path: &Path) -> Result<PathBuf> {
        self.report_progress(
            "database",
            &format!("Preparing database for {fasta_path:?}"),
        );

        // Find FAtoGDB binary
        let fatogdb = self.find_utility("FAtoGDB")?;

        // The output will be .1gdb not .gdb
        let gdb_path = fasta_path.with_extension("1gdb");

        // Check if already exists and is up-to-date
        if gdb_path.exists() {
            let fasta_time = std::fs::metadata(fasta_path)?.modified()?;
            let gdb_time = std::fs::metadata(&gdb_path)?.modified()?;

            if gdb_time >= fasta_time {
                self.report_progress("database", "Using existing database (up-to-date)");
                return Ok(gdb_path);
            }
        }

        // Run FAtoGDB
        self.report_progress("database", "Converting FASTA to GDB format");
        let output = Command::new(&fatogdb).arg(fasta_path).output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(FastGAError::Other(format!("FAtoGDB failed: {stderr}")));
        }

        if !gdb_path.exists() {
            return Err(FastGAError::Other("GDB file was not created".to_string()));
        }

        self.report_progress("database", &format!("Database created: {gdb_path:?}"));
        Ok(gdb_path)
    }

    /// Step 3: Create k-mer index
    pub fn create_index(&self, gdb_path: &Path) -> Result<PathBuf> {
        self.report_progress("index", &format!("Creating index for {gdb_path:?}"));

        let gixmake = self.find_utility("GIXmake")?;
        let gix_path = gdb_path.with_extension("gix");

        // Check if already exists
        if gix_path.exists() {
            self.report_progress("index", "Using existing index");
            return Ok(gix_path);
        }

        // Run GIXmake
        self.report_progress("index", "Building k-mer index");
        let output = Command::new(&gixmake)
            .arg(format!("-T{}", self.config.num_threads))
            // Note: GIXmake doesn't use -f in newer versions (uses -k for k-mer size with default 40)
            .arg(gdb_path)
            .output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(FastGAError::Other(format!("GIXmake failed: {stderr}")));
        }

        self.report_progress("index", &format!("Index created: {gix_path:?}"));
        Ok(gix_path)
    }

    /// Step 4: Run alignment on prepared databases
    pub fn align_databases(&self, query_db: &Path, target_db: &Path) -> Result<String> {
        self.report_progress("align", "Running alignment");

        let fastga = self.find_utility("FastGA")?;
        let alnto_paf = self.find_utility("ALNtoPAF")?;

        // FastGA expects the original FASTA paths, not the GDB paths
        // Extract the original paths by removing the .1gdb extension
        let query_fasta = query_db.with_extension("");
        let target_fasta = target_db.with_extension("");

        // Determine working directory for temp .1aln file
        let working_dir = query_fasta
            .parent()
            .ok_or_else(|| FastGAError::Other("Cannot determine parent directory".to_string()))?;
        let temp_aln = working_dir.join(format!("_tmp_{}.1aln", std::process::id()));

        // Step 1: Run FastGA to create .1aln file
        let mut cmd = Command::new(&fastga);
        cmd.arg(format!(
            "-1:{}",
            temp_aln.file_name().unwrap().to_string_lossy()
        ))
        .arg(format!("-T{}", self.config.num_threads))
        .arg(format!("-l{}", self.config.min_alignment_length));

        if let Some(identity) = self.config.min_identity {
            cmd.arg(format!("-i{identity:.2}"));
        }

        // Use FASTA file names only (FastGA will find .1gdb and .gix in same directory)
        cmd.arg(query_fasta.file_name().unwrap())
            .arg(target_fasta.file_name().unwrap())
            .current_dir(working_dir); // Set working directory

        self.report_progress("align", &format!("Running command: {cmd:?}"));

        let output = cmd.output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let _ = std::fs::remove_file(&temp_aln); // Clean up
            return Err(FastGAError::FastGAExecutionFailed(stderr.to_string()));
        }

        // Step 2: Convert .1aln to PAF using ALNtoPAF
        self.report_progress("align", "Converting .1aln to PAF");
        let paf_cmd = Command::new(&alnto_paf)
            .arg("-x") // Extended CIGAR with X/= operators
            .arg(format!("-T{}", self.config.num_threads))
            .arg(temp_aln.file_name().unwrap())
            .current_dir(working_dir)
            .output()?;

        // Clean up temp .1aln file
        let _ = std::fs::remove_file(&temp_aln);

        if !paf_cmd.status.success() {
            let stderr = String::from_utf8_lossy(&paf_cmd.stderr);
            return Err(FastGAError::Other(format!("ALNtoPAF failed: {stderr}")));
        }

        let paf_output = String::from_utf8_lossy(&paf_cmd.stdout).to_string();
        self.report_progress(
            "align",
            &format!("Alignment complete: {} bytes output", paf_output.len()),
        );

        Ok(paf_output)
    }

    /// Find a utility binary in the build directory
    fn find_utility(&self, name: &str) -> Result<PathBuf> {
        // Try multiple locations
        let locations = vec![
            // In OUT_DIR
            std::env::var("OUT_DIR")
                .ok()
                .map(|d| PathBuf::from(d).join(name)),
            // In current exe's build directory
            std::env::current_exe()
                .ok()
                .and_then(|exe| {
                    exe.parent()
                        .and_then(|p| p.parent())
                        .map(|p| p.join("build"))
                })
                .and_then(|build_dir| {
                    std::fs::read_dir(&build_dir).ok().and_then(|entries| {
                        for entry in entries.flatten() {
                            let name_str = entry.file_name();
                            if name_str.to_string_lossy().starts_with("fastga-rs-") {
                                let util_path = entry.path().join("out").join(name);
                                if util_path.exists() {
                                    return Some(util_path);
                                }
                            }
                        }
                        None
                    })
                }),
            // In target/debug/build
            Some(PathBuf::from("target/debug/build")).and_then(|build_dir| {
                std::fs::read_dir(&build_dir).ok().and_then(|entries| {
                    for entry in entries.flatten() {
                        let name_str = entry.file_name();
                        if name_str.to_string_lossy().starts_with("fastga-rs-") {
                            let util_path = entry.path().join("out").join(name);
                            if util_path.exists() {
                                return Some(util_path);
                            }
                        }
                    }
                    None
                })
            }),
        ];

        for location in locations.into_iter().flatten() {
            if location.exists() {
                self.report_progress("find", &format!("Found {name} at {location:?}"));
                return Ok(location);
            }
        }

        Err(FastGAError::Other(format!("{name} utility not found")))
    }

    /// Full pipeline: validate, prepare, index, and align
    pub fn run_full_pipeline(&self, query: &Path, target: &Path) -> Result<String> {
        // Step 1: Validate
        self.validate_inputs(query, target)?;

        // Step 2: Prepare databases
        let query_db = self.prepare_database(query)?;
        let target_db = self.prepare_database(target)?;

        // Step 3: Create indices
        self.create_index(&query_db)?;
        self.create_index(&target_db)?;

        // Step 4: Align
        self.align_databases(&query_db, &target_db)
    }
}
