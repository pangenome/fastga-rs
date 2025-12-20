//! Rust reimplementation of FastGA's orchestration logic
//!
//! This module orchestrates FastGA utilities (FAtoGDB, GIXmake, FastGA)
//! via system calls instead of FFI to avoid hanging issues.

use crate::binary_finder::find_binary;
use crate::error::{FastGAError, Result};
use std::os::raw::{c_char, c_int};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

/// Global counter for generating unique temp file names
static TEMP_FILE_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Generate a unique temp file name combining PID and atomic counter
fn unique_temp_name() -> String {
    let count = TEMP_FILE_COUNTER.fetch_add(1, Ordering::SeqCst);
    format!("_tmp_{}_{}", std::process::id(), count)
}

// FFI bindings only for the core alignment function
#[link(name = "fastga_main", kind = "static")]
extern "C" {
    #[allow(dead_code)]
    fn fastga_main(argc: c_int, argv: *const *const c_char) -> c_int;
}

/// Pure Rust implementation of FastGA's orchestration
pub struct FastGAOrchestrator {
    pub num_threads: i32,
    pub min_length: i32,
    pub min_identity: f64,
    pub kmer_freq: i32,
    pub temp_dir: String,
}

impl Default for FastGAOrchestrator {
    fn default() -> Self {
        Self::new()
    }
}

impl FastGAOrchestrator {
    pub fn new() -> Self {
        Self {
            num_threads: 8,
            min_length: 100,
            min_identity: 0.7,
            kmer_freq: 10,
            temp_dir: std::env::var("TMPDIR").unwrap_or_else(|_| ".".to_string()),
        }
    }

    /// Main alignment function - orchestrates the entire process
    pub fn align(&self, query_path: &Path, target_path: &Path) -> Result<Vec<u8>> {
        eprintln!("[FastGA] Starting alignment process...");
        eprintln!("[FastGA] Query: {query_path:?}");
        eprintln!("[FastGA] Target: {target_path:?}");
        // Step 1: Convert FASTA files to GDB format if needed
        eprintln!("[FastGA] Step 1: Converting FASTA to GDB format...");
        let _query_gdb = self.prepare_gdb(query_path)?;
        eprintln!("[FastGA] Query GDB created");
        let _target_gdb = self.prepare_gdb(target_path)?;
        eprintln!("[FastGA] Target GDB created");

        // Step 2: Create k-mer indices if needed
        eprintln!("[FastGA] Step 2: Creating k-mer indices...");
        self.create_index(&_query_gdb, self.kmer_freq)?;
        eprintln!("[FastGA] Query index created");
        self.create_index(&_target_gdb, self.kmer_freq)?;
        eprintln!("[FastGA] Target index created");

        // Step 3: Perform the actual alignment (pass FASTA paths, FastGA will use the .gdb/.gix)
        eprintln!("[FastGA] Step 3: Running alignment...");
        let paf_output = self.run_alignment(query_path, target_path)?;
        eprintln!(
            "[FastGA] Alignment complete, output size: {} bytes",
            paf_output.len()
        );

        Ok(paf_output)
    }

    /// Align two genomes assuming GDB/GIX files already exist
    /// This is more efficient when indices are pre-built
    pub fn align_with_existing_indices(
        &self,
        query_path: &Path,
        target_path: &Path,
    ) -> Result<Vec<u8>> {
        self.run_alignment(query_path, target_path)
    }

    /// Align and return .1aln file path directly (NO PAF conversion)
    /// This is the native FastGA output format - most efficient for chaining operations
    pub fn align_to_1aln(&self, query_path: &Path, target_path: &Path) -> Result<String> {
        eprintln!(
            "[FastGA] run_alignment: Aligning {} vs {}",
            query_path.display(),
            target_path.display()
        );

        // Call FastGA binary via system call
        let fastga_bin = find_binary("FastGA")?;

        // Set working directory to where the input files are and use relative paths
        let working_dir = query_path
            .parent()
            .ok_or_else(|| FastGAError::Other("Cannot determine parent directory".to_string()))?;

        let query_filename = query_path
            .file_name()
            .ok_or_else(|| FastGAError::Other("Cannot determine query filename".to_string()))?;
        let target_filename = target_path
            .file_name()
            .ok_or_else(|| FastGAError::Other("Cannot determine target filename".to_string()))?;

        // Create temporary .1aln file with unique name (PID + counter)
        let temp_aln = working_dir.join(format!("{}.1aln", unique_temp_name()));
        let temp_aln_rel = temp_aln.file_name().unwrap();

        // Build command first so we can log actual args
        let mut cmd = std::process::Command::new(&fastga_bin);
        // Set TMPDIR for ONEcode's internal temp files (OneTextSchema, etc.)
        cmd.env("TMPDIR", &self.temp_dir);
        cmd.arg(format!("-1:{}", temp_aln_rel.to_string_lossy()))
            .arg(format!("-T{}", self.num_threads));

        // Add -f for adaptive seed cutoff (FastGA default is 10)
        if self.kmer_freq != 10 {
            cmd.arg(format!("-f{}", self.kmer_freq));
        }

        // Only add -l if it's > 0 (FastGA's default is 0 anyway)
        if self.min_length > 0 {
            cmd.arg(format!("-l{}", self.min_length));
        }

        // Only add -i if it's > 0.0 (let FastGA use its default otherwise)
        if self.min_identity > 0.0 {
            cmd.arg(format!("-i{:.2}", self.min_identity));
        }

        cmd.arg(query_filename)
            .arg(target_filename)
            .current_dir(working_dir);

        // Log the actual command being executed
        let args_str: Vec<String> = cmd
            .get_args()
            .map(|s| s.to_string_lossy().to_string())
            .collect();
        eprintln!(
            "[FastGA] Calling FastGA: {} {} (in dir: {}, TMPDIR={})",
            fastga_bin.display(),
            args_str.join(" "),
            working_dir.display(),
            self.temp_dir
        );

        let output = cmd
            .output()
            .map_err(|e| FastGAError::Other(format!("Failed to run FastGA: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stdout = String::from_utf8_lossy(&output.stdout);
            // Clean up temp file
            let _ = std::fs::remove_file(&temp_aln);
            return Err(FastGAError::Other(format!(
                "FastGA failed with code {:?}\nstdout: {}\nstderr: {}",
                output.status.code(),
                stdout,
                stderr
            )));
        }

        eprintln!(
            "[FastGA] FastGA completed, returning .1aln file: {}",
            temp_aln.display()
        );
        Ok(temp_aln.to_string_lossy().to_string())
    }

    /// Convert FASTA to GDB format using FAtoGDB
    pub fn prepare_gdb(&self, fasta_path: &Path) -> Result<String> {
        eprintln!("[FastGA] prepare_gdb: Converting {fasta_path:?} to GDB");
        let gdb_base = fasta_path.with_extension("");
        let gdb_path = format!("{}.gdb", gdb_base.display());
        let gdb_path_1 = format!("{}.1gdb", gdb_base.display());

        // Check if GDB already exists and is up-to-date
        if Path::new(&gdb_path).exists() || Path::new(&gdb_path_1).exists() {
            if let Ok(fasta_meta) = std::fs::metadata(fasta_path) {
                let gdb_exists = Path::new(&gdb_path).exists();
                let gdb_1_exists = Path::new(&gdb_path_1).exists();

                if gdb_exists || gdb_1_exists {
                    let check_path = if gdb_1_exists { &gdb_path_1 } else { &gdb_path };
                    if let Ok(gdb_meta) = std::fs::metadata(check_path) {
                        if let (Ok(fasta_time), Ok(gdb_time)) =
                            (fasta_meta.modified(), gdb_meta.modified())
                        {
                            if gdb_time >= fasta_time {
                                eprintln!("[FastGA] GDB already exists and is up-to-date");
                                return Ok(gdb_base.to_string_lossy().to_string());
                            }
                        }
                    }
                }
            }
        }

        // Call FAtoGDB via system call
        let fatogdb_bin = find_binary("FAtoGDB")?;

        eprintln!(
            "[FastGA] Calling FAtoGDB: {} {} (TMPDIR={})",
            fatogdb_bin.display(),
            fasta_path.display(),
            self.temp_dir
        );

        let output = std::process::Command::new(&fatogdb_bin)
            .env("TMPDIR", &self.temp_dir)
            .arg(fasta_path)
            .output()
            .map_err(|e| FastGAError::Other(format!("Failed to run FAtoGDB: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stdout = String::from_utf8_lossy(&output.stdout);
            return Err(FastGAError::Other(format!(
                "FAtoGDB failed with code {:?}\nstdout: {}\nstderr: {}",
                output.status.code(),
                stdout,
                stderr
            )));
        }

        eprintln!("[FastGA] FAtoGDB completed successfully");
        Ok(gdb_base.to_string_lossy().to_string())
    }

    /// Create k-mer index using GIXmake
    pub fn create_index(&self, gdb_path: &str, _freq: i32) -> Result<()> {
        eprintln!("[FastGA] create_index: Creating index for {gdb_path}");
        let gix_path = format!("{gdb_path}.gix");

        // Check if index already exists
        if Path::new(&gix_path).exists() {
            eprintln!("[FastGA] Index already exists: {gix_path}");
            return Ok(());
        }

        // Call GIXmake via system call
        // Note: GIXmake doesn't use -f for frequency, it uses -k for k-mer size (default 40)
        let gixmake_bin = find_binary("GIXmake")?;

        eprintln!(
            "[FastGA] Calling GIXmake: {} -T{} -P{} {}",
            gixmake_bin.display(),
            self.num_threads,
            self.temp_dir,
            gdb_path
        );

        let output = std::process::Command::new(&gixmake_bin)
            .env("TMPDIR", &self.temp_dir)
            .arg(format!("-T{}", self.num_threads))
            .arg(format!("-P{}", self.temp_dir))
            .arg(gdb_path)
            .output()
            .map_err(|e| FastGAError::Other(format!("Failed to run GIXmake: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stdout = String::from_utf8_lossy(&output.stdout);
            return Err(FastGAError::Other(format!(
                "GIXmake failed with code {:?}\nstdout: {}\nstderr: {}",
                output.status.code(),
                stdout,
                stderr
            )));
        }

        eprintln!("[FastGA] GIXmake completed successfully");
        Ok(())
    }

    /// Run the actual alignment algorithm using FastGA binary
    fn run_alignment(&self, query_path: &Path, target_path: &Path) -> Result<Vec<u8>> {
        eprintln!(
            "[FastGA] run_alignment: Aligning {} vs {}",
            query_path.display(),
            target_path.display()
        );

        // Call FastGA binary via system call
        let fastga_bin = find_binary("FastGA")?;
        let alnto_paf_bin = find_binary("ALNtoPAF")?;

        // Set working directory to where the input files are and use relative paths
        let working_dir = query_path
            .parent()
            .ok_or_else(|| FastGAError::Other("Cannot determine parent directory".to_string()))?;

        let query_filename = query_path
            .file_name()
            .ok_or_else(|| FastGAError::Other("Cannot determine query filename".to_string()))?;
        let target_filename = target_path
            .file_name()
            .ok_or_else(|| FastGAError::Other("Cannot determine target filename".to_string()))?;

        // Create temporary .1aln file with unique name (PID + counter)
        let temp_aln = working_dir.join(format!("{}.1aln", unique_temp_name()));
        let temp_aln_rel = temp_aln.file_name().unwrap();

        // Build command first so we can log actual args
        let mut cmd = std::process::Command::new(&fastga_bin);
        // Set TMPDIR for ONEcode's internal temp files (OneTextSchema, etc.)
        cmd.env("TMPDIR", &self.temp_dir);
        cmd.arg(format!("-1:{}", temp_aln_rel.to_string_lossy()))
            .arg(format!("-T{}", self.num_threads));

        // Add -f for adaptive seed cutoff (FastGA default is 10)
        if self.kmer_freq != 10 {
            cmd.arg(format!("-f{}", self.kmer_freq));
        }

        // Only add -l if it's > 0 (FastGA's default is 0 anyway)
        if self.min_length > 0 {
            cmd.arg(format!("-l{}", self.min_length));
        }

        // Only add -i if it's > 0.0 (let FastGA use its default otherwise)
        if self.min_identity > 0.0 {
            cmd.arg(format!("-i{:.2}", self.min_identity));
        }

        cmd.arg(query_filename)
            .arg(target_filename)
            .current_dir(working_dir);

        // Log the actual command being executed
        let args_str: Vec<String> = cmd
            .get_args()
            .map(|s| s.to_string_lossy().to_string())
            .collect();
        eprintln!(
            "[FastGA] Calling FastGA: {} {} (in dir: {}, TMPDIR={})",
            fastga_bin.display(),
            args_str.join(" "),
            working_dir.display(),
            self.temp_dir
        );

        let output = cmd
            .output()
            .map_err(|e| FastGAError::Other(format!("Failed to run FastGA: {e}")))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stdout = String::from_utf8_lossy(&output.stdout);
            // Clean up temp file
            let _ = std::fs::remove_file(&temp_aln);
            return Err(FastGAError::Other(format!(
                "FastGA failed with code {:?}\nstdout: {}\nstderr: {}",
                output.status.code(),
                stdout,
                stderr
            )));
        }

        eprintln!("[FastGA] FastGA completed, converting .1aln to PAF");

        // Now convert .1aln to PAF using ALNtoPAF
        // Using -x flag for extended CIGAR with X/= operators (works with FastGA 1789f15)
        let paf_output = std::process::Command::new(&alnto_paf_bin)
            .env("TMPDIR", &self.temp_dir)
            .arg("-x")
            .arg(temp_aln_rel)
            .current_dir(working_dir)
            .output()
            .map_err(|e| FastGAError::Other(format!("Failed to run ALNtoPAF: {e}")))?;

        // Clean up temp .1aln file
        let _ = std::fs::remove_file(&temp_aln);

        if !paf_output.status.success() {
            let stderr = String::from_utf8_lossy(&paf_output.stderr);
            let stdout = String::from_utf8_lossy(&paf_output.stdout);
            return Err(FastGAError::Other(format!(
                "ALNtoPAF failed with code {:?}\nstdout: {}\nstderr: {}",
                paf_output.status.code(),
                stdout,
                stderr
            )));
        }

        eprintln!(
            "[FastGA] Conversion completed successfully, output size: {} bytes",
            paf_output.stdout.len()
        );
        Ok(paf_output.stdout)
    }

    /// Clean up temporary files
    #[allow(dead_code)]
    fn cleanup(&self, _gdb_path: &str) -> Result<()> {
        // Optionally remove GDB and GIX files
        // For now, leave them as cache
        Ok(())
    }
}

/// Direct FFI-based alignment that bypasses FastGA's main entirely
pub fn align_direct(query: &Path, target: &Path, config: crate::Config) -> Result<Vec<u8>> {
    let orchestrator = FastGAOrchestrator {
        num_threads: config.num_threads as i32,
        min_length: config.min_alignment_length as i32,
        min_identity: config.min_identity.unwrap_or(0.7),
        kmer_freq: 10,
        temp_dir: std::env::var("TMPDIR").unwrap_or_else(|_| ".".to_string()),
    };

    orchestrator.align(query, target)
}
