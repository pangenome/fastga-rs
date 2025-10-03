// Native .1aln reader - direct reading without PAF conversion
//
// This module provides efficient reading of .1aln alignment files directly
// using the ONEcode library, avoiding the overhead of converting to PAF format.

use anyhow::Result;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void};

/// Alignment record from .1aln file
#[repr(C)]
#[derive(Debug, Clone)]
pub struct AlnRecord {
    pub query_id: i64,
    pub query_start: i64,
    pub query_end: i64,
    pub target_id: i64,
    pub target_start: i64,
    pub target_end: i64,
    pub reverse: c_int,
    pub diffs: c_int,
    pub query_len: i64,
    pub target_len: i64,
}

extern "C" {
    fn aln_open(path: *const c_char, num_alignments: *mut i64) -> *mut c_void;
    fn aln_read_record(handle: *mut c_void, record: *mut AlnRecord) -> c_int;
    fn aln_get_seq_name(handle: *mut c_void, seq_id: i64, which_db: c_int) -> *const c_char;
    fn aln_close(handle: *mut c_void);

    // Writer functions
    fn aln_create(path: *const c_char, gdb1_path: *const c_char, gdb2_path: *const c_char) -> *mut c_void;
    fn aln_write_record(handle: *mut c_void, record: *const AlnRecord) -> c_int;
    fn aln_close_writer(handle: *mut c_void);
}

/// Reader for .1aln alignment files
pub struct AlnReader {
    handle: *mut c_void,
    num_alignments: i64,
}

impl AlnReader {
    /// Open a .1aln file for reading
    ///
    /// # Requirements
    /// - The .1aln file must exist
    /// - Matching .1gdb file(s) must be present
    /// - Hidden .bps file(s) must be present
    ///
    /// # Example
    /// ```no_run
    /// use fastga_rs::AlnReader;
    ///
    /// let mut reader = AlnReader::open("alignments.1aln")?;
    /// println!("Total alignments: {}", reader.num_alignments());
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn open(path: &str) -> Result<Self> {
        let c_path = CString::new(path)?;
        let mut num_alignments: i64 = 0;

        let handle = unsafe { aln_open(c_path.as_ptr(), &mut num_alignments) };

        if handle.is_null() {
            anyhow::bail!("Failed to open .1aln file: {}. Check that matching .1gdb and .bps files exist.", path);
        }

        Ok(AlnReader {
            handle,
            num_alignments,
        })
    }

    /// Get total number of alignments in the file
    pub fn num_alignments(&self) -> i64 {
        self.num_alignments
    }

    /// Read the next alignment record
    ///
    /// Returns `Ok(Some(record))` if a record was read, `Ok(None)` at EOF
    ///
    /// # Example
    /// ```no_run
    /// use fastga_rs::AlnReader;
    ///
    /// let mut reader = AlnReader::open("alignments.1aln")?;
    /// while let Some(rec) = reader.read_record()? {
    ///     println!("Query: {}..{}", rec.query_start, rec.query_end);
    /// }
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn read_record(&mut self) -> Result<Option<AlnRecord>> {
        let mut record = AlnRecord {
            query_id: 0,
            query_start: 0,
            query_end: 0,
            target_id: 0,
            target_start: 0,
            target_end: 0,
            reverse: 0,
            diffs: 0,
            query_len: 0,
            target_len: 0,
        };

        let result = unsafe { aln_read_record(self.handle, &mut record) };

        if result < 0 {
            Ok(None) // EOF
        } else {
            Ok(Some(record))
        }
    }

    /// Get sequence name by scaffold ID
    ///
    /// # Arguments
    /// - `seq_id`: Scaffold ID from AlnRecord
    /// - `which_db`: 0 for query database, 1 for target database
    ///
    /// # Example
    /// ```no_run
    /// use fastga_rs::AlnReader;
    ///
    /// let mut reader = AlnReader::open("alignments.1aln")?;
    /// if let Some(rec) = reader.read_record()? {
    ///     let qname = reader.get_seq_name(rec.query_id, 0)?;
    ///     let tname = reader.get_seq_name(rec.target_id, 1)?;
    ///     println!("{} -> {}", qname, tname);
    /// }
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn get_seq_name(&self, seq_id: i64, which_db: i32) -> Result<String> {
        let c_str = unsafe { aln_get_seq_name(self.handle, seq_id, which_db) };

        if c_str.is_null() {
            anyhow::bail!("Failed to get sequence name for ID {}", seq_id);
        }

        let name = unsafe { CStr::from_ptr(c_str) }
            .to_str()?
            .to_string();

        Ok(name)
    }
}

impl Drop for AlnReader {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { aln_close(self.handle) };
        }
    }
}

// SAFETY: AlnReader handle is only used from one thread (not shared)
// The C library doesn't maintain global state
unsafe impl Send for AlnReader {}

/// Writer for .1aln alignment files
pub struct AlnWriter {
    handle: *mut c_void,
}

impl AlnWriter {
    /// Create a new .1aln file for writing
    ///
    /// # Arguments
    /// - `path`: Output .1aln file path
    /// - `gdb1_path`: Path to query .1gdb file (sequence metadata)
    /// - `gdb2_path`: Path to target .1gdb file (sequence metadata)
    ///
    /// # Example
    /// ```no_run
    /// use fastga_rs::AlnWriter;
    ///
    /// let mut writer = AlnWriter::create(
    ///     "output.1aln",
    ///     "sequences.1gdb",
    ///     "sequences.1gdb"
    /// )?;
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn create(path: &str, gdb1_path: &str, gdb2_path: &str) -> Result<Self> {
        let c_path = CString::new(path)?;
        let c_gdb1 = CString::new(gdb1_path)?;
        let c_gdb2 = CString::new(gdb2_path)?;

        let handle = unsafe {
            aln_create(c_path.as_ptr(), c_gdb1.as_ptr(), c_gdb2.as_ptr())
        };

        if handle.is_null() {
            anyhow::bail!("Failed to create .1aln file: {}", path);
        }

        Ok(AlnWriter { handle })
    }

    /// Write an alignment record
    ///
    /// # Example
    /// ```no_run
    /// use fastga_rs::{AlnWriter, AlnRecord};
    ///
    /// let mut writer = AlnWriter::create("out.1aln", "db.1gdb", "db.1gdb")?;
    ///
    /// let rec = AlnRecord {
    ///     query_id: 0,
    ///     query_start: 0,
    ///     query_end: 1000,
    ///     target_id: 1,
    ///     target_start: 500,
    ///     target_end: 1500,
    ///     reverse: 0,
    ///     diffs: 10,
    ///     query_len: 5000,
    ///     target_len: 6000,
    /// };
    ///
    /// writer.write_record(&rec)?;
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &AlnRecord) -> Result<()> {
        let result = unsafe { aln_write_record(self.handle, record) };

        if result < 0 {
            anyhow::bail!("Failed to write alignment record");
        }

        Ok(())
    }
}

impl Drop for AlnWriter {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { aln_close_writer(self.handle) };
        }
    }
}

// SAFETY: AlnWriter handle is only used from one thread
unsafe impl Send for AlnWriter {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore] // Only run when .1aln test file is available
    fn test_aln_reader_basic() {
        let mut reader = AlnReader::open("test.1aln").expect("Failed to open test file");

        println!("Total alignments: {}", reader.num_alignments());

        let mut count = 0;
        while let Some(record) = reader.read_record().expect("Failed to read record") {
            let query_name = reader.get_seq_name(record.query_id, 0).unwrap();
            let target_name = reader.get_seq_name(record.target_id, 1).unwrap();

            println!("Alignment {}: {} -> {}", count, query_name, target_name);
            count += 1;

            if count >= 5 {
                break;
            }
        }

        assert!(count > 0, "Should have read at least one record");
    }

    #[test]
    #[ignore] // Only run when .1aln test file is available
    fn test_native_vs_alnTopaf_coordinates() {
        // This test validates that native .1aln reading produces identical
        // coordinates to ALNtoPAF for filtering purposes

        let mut reader = AlnReader::open("test.1aln").expect("Failed to open test file");

        // Read first record
        let rec = reader.read_record().expect("Failed to read first record")
            .expect("No records in file");

        // Expected coordinates from ALNtoPAF (first record):
        // SGDref#1#chrI	230218	0	12550	-	DBVPG6765#1#chrVIII	533397	520779	533304

        // Validate query coordinates
        assert_eq!(rec.query_start, 0, "Query start should be 0");
        assert_eq!(rec.query_end, 12550, "Query end should be 12550");
        assert_eq!(rec.query_len, 230218, "Query length should be 230218");

        // Validate target coordinates (reverse strand, so adjusted)
        assert_eq!(rec.target_start, 520779, "Target start should be 520779");
        assert_eq!(rec.target_end, 533304, "Target end should be 533304");
        assert_eq!(rec.target_len, 533397, "Target length should be 533397");

        // Validate reverse flag
        assert_eq!(rec.reverse, 1, "Should be reverse strand");

        // Validate sequence names
        let qname = reader.get_seq_name(rec.query_id, 0).unwrap();
        let tname = reader.get_seq_name(rec.target_id, 1).unwrap();

        assert_eq!(qname, "SGDref#1#chrI", "Query name mismatch");
        assert_eq!(tname, "DBVPG6765#1#chrVIII", "Target name mismatch");

        println!("✓ First record coordinates match ALNtoPAF exactly");
    }

    #[test]
    #[ignore] // Only run when .1aln test file is available
    fn test_coordinate_properties() {
        // Test that all coordinates satisfy required properties
        let mut reader = AlnReader::open("test.1aln").expect("Failed to open test file");

        let mut count = 0;
        while let Some(rec) = reader.read_record().expect("Failed to read record") {
            // Property 1: Coordinates are non-negative
            assert!(rec.query_start >= 0, "Negative query start at record {}", count);
            assert!(rec.query_end >= 0, "Negative query end at record {}", count);
            assert!(rec.target_start >= 0, "Negative target start at record {}", count);
            assert!(rec.target_end >= 0, "Negative target end at record {}", count);

            // Property 2: Start < End
            assert!(rec.query_start < rec.query_end,
                    "Invalid query interval at record {}: {}..{}",
                    count, rec.query_start, rec.query_end);
            assert!(rec.target_start < rec.target_end,
                    "Invalid target interval at record {}: {}..{}",
                    count, rec.target_start, rec.target_end);

            // Property 3: Coordinates within sequence bounds
            assert!(rec.query_end <= rec.query_len,
                    "Query end {} > length {} at record {}",
                    rec.query_end, rec.query_len, count);
            assert!(rec.target_end <= rec.target_len,
                    "Target end {} > length {} at record {}",
                    rec.target_end, rec.target_len, count);

            count += 1;

            if count >= 1000 {
                break; // Test first 1000 records
            }
        }

        assert!(count > 0, "Should have read at least one record");
        println!("✓ Validated properties for {} records", count);
    }
}
