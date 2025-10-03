/// ONElib FFI bindings for reading/writing .1aln files
///
/// This provides safe Rust wrappers around the ONElib C library
/// for working with .1aln alignment files.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};
use std::path::Path;
use std::ptr;
use anyhow::{Result, bail, Context};

use crate::alignment::Alignment;

// ONElib C types and functions
#[repr(C)]
struct OneFile {
    _private: [u8; 0], // Opaque type
}

#[repr(C)]
struct OneSchema {
    _private: [u8; 0], // Opaque type
}

#[repr(C)]
union OneField {
    i: i64,
    r: f64,
    c: c_char,
    len: i64,
}

extern "C" {
    // Schema creation
    fn make_Aln_Schema() -> *mut OneSchema;
    fn oneSchemaDestroy(schema: *mut OneSchema);

    // File operations
    fn oneFileOpenRead(
        path: *const c_char,
        schema: *mut OneSchema,
        file_type: *const c_char,
        nthreads: c_int,
    ) -> *mut OneFile;

    fn oneFileOpenWriteNew(
        path: *const c_char,
        schema: *mut OneSchema,
        file_type: *const c_char,
        is_binary: bool,
        nthreads: c_int,
    ) -> *mut OneFile;

    fn oneFileClose(file: *mut OneFile);

    // Reading
    fn oneReadLine(file: *mut OneFile) -> c_char;

    // Field accessors (from rust_onelib.c)
    fn one_int(file: *mut OneFile, index: c_int) -> i64;
    fn one_real(file: *mut OneFile, index: c_int) -> f64;
    fn one_char(file: *mut OneFile, index: c_int) -> c_char;
    fn one_line_type(file: *mut OneFile) -> c_char;
    fn one_line_count(file: *mut OneFile) -> i64;
}

/// Reader for .1aln files
pub struct AlnReader {
    file: *mut OneFile,
    schema: *mut OneSchema,
}

impl AlnReader {
    /// Open a .1aln file for reading
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_str = path.as_ref().to_str()
            .context("Invalid path")?;
        let path_cstr = CString::new(path_str)?;

        unsafe {
            let schema = make_Aln_Schema();
            if schema.is_null() {
                bail!("Failed to create .1aln schema");
            }

            let type_cstr = CString::new("aln")?;
            let file = oneFileOpenRead(
                path_cstr.as_ptr(),
                schema,
                type_cstr.as_ptr(),
                1, // single threaded for now
            );

            if file.is_null() {
                oneSchemaDestroy(schema);
                bail!("Failed to open .1aln file: {}", path_str);
            }

            Ok(AlnReader { file, schema })
        }
    }

    /// Read next alignment from the file
    /// Returns None when EOF is reached
    pub fn read_alignment(&mut self) -> Result<Option<Alignment>> {
        unsafe {
            // Skip to next 'A' record (alignment)
            loop {
                let line_type = one_line_type(self.file);

                // If we're not at an 'A' record, read next line
                if line_type != b'A' as c_char {
                    let next = oneReadLine(self.file);
                    if next == 0 {
                        return Ok(None); // EOF
                    }
                    continue;
                }

                // We're at an 'A' record
                // Schema: O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
                // Fields: aread, abpos, aepos, bread, bbpos, bepos

                let a_id = one_int(self.file, 0);
                let a_beg = one_int(self.file, 1) as u64;
                let a_end = one_int(self.file, 2) as u64;
                let b_id = one_int(self.file, 3);
                let b_beg = one_int(self.file, 4) as u64;
                let b_end = one_int(self.file, 5) as u64;

                // Read associated data lines until we hit 'T' (trace) or next 'A'
                let mut matches = 0u64;
                let mut diffs = 0u64;
                let mut is_reverse = false;
                let mut query_len = 0u64;
                let mut target_len = 0u64;

                loop {
                    let next_type = oneReadLine(self.file);

                    if next_type == 0 {
                        break; // EOF
                    }

                    match next_type as u8 as char {
                        'T' => {
                            // Trace record - marks end of this alignment's data
                            // Read next line to position for next alignment
                            oneReadLine(self.file);
                            break;
                        }
                        'M' => matches = one_int(self.file, 0) as u64,
                        'D' => diffs = one_int(self.file, 0) as u64,
                        'R' => is_reverse = true,
                        'L' => {
                            query_len = one_int(self.file, 0) as u64;
                            target_len = one_int(self.file, 1) as u64;
                        }
                        'A' => {
                            // Hit next alignment without seeing 'T'
                            // This is valid - not all alignments have traces
                            break;
                        }
                        _ => {
                            // Skip other records (X, C, Q, etc.)
                            continue;
                        }
                    }
                }

                // Calculate block length and identity
                let block_len = a_end - a_beg;
                let alignment_len = matches + diffs;

                // Identity: matches / total_aligned_bases
                let identity = if alignment_len > 0 {
                    matches as f64 / alignment_len as f64
                } else {
                    0.0
                };

                // Create alignment using sequence IDs
                // Filtering doesn't need actual names, just unique identifiers
                let alignment = Alignment {
                    query_name: format!("{}", a_id),
                    query_len: query_len as usize,
                    query_start: a_beg as usize,
                    query_end: a_end as usize,
                    strand: if is_reverse { '-' } else { '+' },
                    target_name: format!("{}", b_id),
                    target_len: target_len as usize,
                    target_start: b_beg as usize,
                    target_end: b_end as usize,
                    matches: matches as usize,
                    block_len: block_len as usize,
                    mapping_quality: 60, // Default quality
                    cigar: String::new(), // Don't extract CIGAR for filtering
                    tags: Vec::new(),
                    mismatches: diffs as usize, // Differences include substitutions + indels
                    gap_opens: 0, // Not tracked in basic .1aln
                    gap_len: 0,
                };

                return Ok(Some(alignment));
            }
        }
    }

    /// Read all alignments from the file
    pub fn read_all(&mut self) -> Result<Vec<Alignment>> {
        let mut alignments = Vec::new();

        while let Some(alignment) = self.read_alignment()? {
            alignments.push(alignment);
        }

        Ok(alignments)
    }
}

impl Drop for AlnReader {
    fn drop(&mut self) {
        unsafe {
            if !self.file.is_null() {
                oneFileClose(self.file);
            }
            if !self.schema.is_null() {
                oneSchemaDestroy(self.schema);
            }
        }
    }
}

/// Writer for .1aln files
pub struct AlnWriter {
    file: *mut OneFile,
    schema: *mut OneSchema,
}

impl AlnWriter {
    /// Create a new .1aln file for writing
    pub fn create<P: AsRef<Path>>(path: P, binary: bool) -> Result<Self> {
        let path_str = path.as_ref().to_str()
            .context("Invalid path")?;
        let path_cstr = CString::new(path_str)?;

        unsafe {
            let schema = make_Aln_Schema();
            if schema.is_null() {
                bail!("Failed to create .1aln schema");
            }

            let type_cstr = CString::new("aln")?;
            let file = oneFileOpenWriteNew(
                path_cstr.as_ptr(),
                schema,
                type_cstr.as_ptr(),
                binary,
                1, // single threaded
            );

            if file.is_null() {
                oneSchemaDestroy(schema);
                bail!("Failed to create .1aln file: {}", path_str);
            }

            Ok(AlnWriter { file, schema })
        }
    }

    /// Write an alignment to the file
    pub fn write_alignment(&mut self, _aln: &Alignment) -> Result<()> {
        // TODO: Implement writing
        // This is more complex as we need to write multiple line types
        bail!("Writing .1aln not yet implemented")
    }
}

impl Drop for AlnWriter {
    fn drop(&mut self) {
        unsafe {
            if !self.file.is_null() {
                oneFileClose(self.file);
            }
            if !self.schema.is_null() {
                oneSchemaDestroy(self.schema);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore] // Only run with actual .1aln files
    fn test_read_1aln() {
        let mut reader = AlnReader::open("test.1aln").unwrap();
        let alignments = reader.read_all().unwrap();
        assert!(alignments.len() > 0);
    }
}
