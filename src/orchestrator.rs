//! Rust reimplementation of FastGA's orchestration logic
//!
//! This module replaces FastGA's main() with pure Rust code that calls
//! the underlying C functions directly via FFI.

use crate::error::{FastGAError, Result};
use std::ffi::CString;
use std::os::raw::{c_char, c_int};
use std::path::Path;

// FFI bindings to the actual worker functions in FastGA's C code
#[link(name = "fatogdb_main", kind = "static")]
#[link(name = "gixmake_main", kind = "static")]
#[link(name = "fastga_main", kind = "static")]
#[link(name = "fastga_common", kind = "static")]
extern "C" {
    // FAtoGDB functions
    fn fatogdb_main(argc: c_int, argv: *const *const c_char) -> c_int;

    // GIXmake functions
    fn gixmake_main(argc: c_int, argv: *const *const c_char) -> c_int;

    // GDB functions would go here if we needed them directly
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
        let query_gdb = self.prepare_gdb(query_path)?;
        eprintln!("[FastGA] Query GDB created: {query_gdb}");
        let target_gdb = self.prepare_gdb(target_path)?;
        eprintln!("[FastGA] Target GDB created: {target_gdb}");

        // Step 2: Create k-mer indices if needed
        eprintln!("[FastGA] Step 2: Creating k-mer indices...");
        self.create_index(&query_gdb, self.kmer_freq)?;
        eprintln!("[FastGA] Query index created");
        self.create_index(&target_gdb, self.kmer_freq)?;
        eprintln!("[FastGA] Target index created");

        // Step 3: Perform the actual alignment
        eprintln!("[FastGA] Step 3: Running alignment...");
        let paf_output = self.run_alignment(&query_gdb, &target_gdb)?;
        eprintln!("[FastGA] Alignment complete, output size: {} bytes", paf_output.len());

        // Clean up temporary files
        self.cleanup(&query_gdb)?;
        self.cleanup(&target_gdb)?;

        Ok(paf_output)
    }

    /// Convert FASTA to GDB format using FAtoGDB
    fn prepare_gdb(&self, fasta_path: &Path) -> Result<String> {
        eprintln!("[FastGA] prepare_gdb: Converting {fasta_path:?} to GDB");
        let gdb_path = fasta_path.with_extension("gdb");

        // Check if GDB already exists and is up-to-date
        if gdb_path.exists() {
            if let (Ok(fasta_meta), Ok(gdb_meta)) =
                (std::fs::metadata(fasta_path), std::fs::metadata(&gdb_path))
            {
                if let (Ok(fasta_time), Ok(gdb_time)) = (fasta_meta.modified(), gdb_meta.modified())
                {
                    if gdb_time >= fasta_time {
                        return Ok(gdb_path.to_string_lossy().to_string());
                    }
                }
            }
        }

        // Call FAtoGDB directly via FFI
        let args = [CString::new("FAtoGDB").unwrap(),
            CString::new(fasta_path.to_string_lossy().to_string()).unwrap()];

        let argv: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();

        eprintln!("[FastGA] Calling fatogdb_main with args: {:?}", args.iter().map(|s| s.to_str().unwrap()).collect::<Vec<_>>());
        let result = unsafe { fatogdb_main(argv.len() as c_int, argv.as_ptr()) };
        eprintln!("[FastGA] fatogdb_main returned: {result}");

        if result != 0 {
            return Err(FastGAError::Other(format!(
                "FAtoGDB failed with code {result}"
            )));
        }

        Ok(gdb_path.to_string_lossy().to_string())
    }

    /// Create k-mer index using GIXmake
    fn create_index(&self, gdb_path: &str, freq: i32) -> Result<()> {
        eprintln!("[FastGA] create_index: Creating index for {gdb_path}");
        let gix_path = gdb_path.replace(".gdb", ".gix");

        // Check if index already exists
        if Path::new(&gix_path).exists() {
            // TODO: Check if it's up-to-date with the GDB file
            return Ok(());
        }

        // Call GIXmake directly via FFI
        let args = [CString::new("GIXmake").unwrap(),
            CString::new(format!("-T{}", self.num_threads)).unwrap(),
            CString::new(format!("-P{}", self.temp_dir)).unwrap(),
            CString::new(format!("-f{freq}")).unwrap(),
            CString::new(gdb_path).unwrap()];

        let argv: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();

        eprintln!("[FastGA] Calling gixmake_main with args: {:?}", args.iter().map(|s| s.to_str().unwrap()).collect::<Vec<_>>());
        let result = unsafe { gixmake_main(argv.len() as c_int, argv.as_ptr()) };
        eprintln!("[FastGA] gixmake_main returned: {result}");

        if result != 0 {
            return Err(FastGAError::Other(format!(
                "GIXmake failed with code {result}"
            )));
        }

        Ok(())
    }

    /// Run the actual alignment algorithm
    fn run_alignment(&self, query_gdb: &str, target_gdb: &str) -> Result<Vec<u8>> {
        eprintln!("[FastGA] run_alignment: Aligning {query_gdb} vs {target_gdb}");
        // This is where we'd call the core alignment functions from FastGA
        // For now, we'll still use the wrapped main function, but ideally
        // we'd extract and call the actual alignment logic directly

        use std::os::unix::io::FromRawFd;

        // Create pipe to capture output
        let (read_fd, write_fd) = nix::unistd::pipe()
            .map_err(|e| FastGAError::Other(format!("Failed to create pipe: {e}")))?;

        // Save original stdout
        let original_stdout = unsafe { libc::dup(1) };

        // Redirect stdout to pipe
        unsafe {
            libc::dup2(write_fd, 1);
            libc::close(write_fd);
        }

        // Build arguments for the core alignment
        // We're calling the main function but with GDB files already prepared
        let args = [CString::new("FastGA").unwrap(),
            CString::new("-pafx").unwrap(),
            CString::new(format!("-T{}", self.num_threads)).unwrap(),
            CString::new(format!("-l{}", self.min_length)).unwrap(),
            CString::new(format!("-i{:.2}", self.min_identity)).unwrap(),
            CString::new(query_gdb).unwrap(),
            CString::new(target_gdb).unwrap()];

        let argv: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();

        // For now, still calling fastga_main, but we should extract the actual
        // alignment logic and call it directly
        extern "C" {
            fn fastga_main(argc: c_int, argv: *const *const c_char) -> c_int;
        }

        eprintln!("[FastGA] Calling fastga_main with args: {:?}", args.iter().map(|s| s.to_str().unwrap()).collect::<Vec<_>>());
        let result = unsafe { fastga_main(argv.len() as c_int, argv.as_ptr()) };
        eprintln!("[FastGA] fastga_main returned: {result}");

        // Restore stdout
        unsafe {
            libc::dup2(original_stdout, 1);
            libc::close(original_stdout);
        }

        if result != 0 {
            return Err(FastGAError::Other(format!(
                "Alignment failed with code {result}"
            )));
        }

        // Read output from pipe
        let mut output = Vec::new();
        let mut pipe_file = unsafe { std::fs::File::from_raw_fd(read_fd) };
        std::io::Read::read_to_end(&mut pipe_file, &mut output)
            .map_err(FastGAError::IoError)?;

        Ok(output)
    }

    /// Clean up temporary files
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
