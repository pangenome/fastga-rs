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
            // Read lines until we find an 'A' (alignment) record
            loop {
                let line_type = oneReadLine(self.file);

                if line_type == 0 {
                    // EOF
                    return Ok(None);
                }

                if line_type as u8 as char == 'A' {
                    // Found alignment record
                    // Schema: O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
                    //         positions: a_id, a_beg, a_end, b_id, b_beg, b_end

                    let a_id = one_int(self.file, 0) as usize;
                    let a_beg = one_int(self.file, 1) as u64;
                    let a_end = one_int(self.file, 2) as u64;
                    let b_id = one_int(self.file, 3) as usize;
                    let b_beg = one_int(self.file, 4) as u64;
                    let b_end = one_int(self.file, 5) as u64;

                    // Read associated data lines (M, D, R, etc.)
                    let mut matches = 0;
                    let mut diffs = 0;
                    let mut is_reverse = false;
                    let mut query_len = 0;
                    let mut target_len = 0;

                    loop {
                        let next_type = oneReadLine(self.file);

                        if next_type == 0 || next_type as u8 as char == 'A' {
                            // EOF or next alignment - we're done with this one
                            // Note: We need to "put back" this line somehow
                            // For now, just process what we have
                            break;
                        }

                        match next_type as u8 as char {
                            'M' => matches = one_int(self.file, 0) as u64,
                            'D' => diffs = one_int(self.file, 0) as u64,
                            'R' => is_reverse = true,
                            'L' => {
                                query_len = one_int(self.file, 0) as u64;
                                target_len = one_int(self.file, 1) as u64;
                            }
                            'T' | 'X' | 'C' | 'Q' => {
                                // Skip trace points, diff counts, CIGAR, quality
                                // We don't need these for filtering
                                continue;
                            }
                            _ => continue,
                        }
                    }

                    // Calculate identity
                    let block_len = a_end - a_beg;
                    let identity = if block_len > 0 {
                        matches as f64 / (matches + diffs) as f64
                    } else {
                        0.0
                    };

                    // Create alignment
                    // Note: We don't have sequence names yet - those need to come from
                    // the sequence section of the file. For now, use IDs as placeholder.
                    let alignment = Alignment {
                        query_name: format!("seq_{}", a_id),
                        query_len,
                        query_start: a_beg,
                        query_end: a_end,
                        strand: if is_reverse { '-' } else { '+' },
                        target_name: format!("seq_{}", b_id),
                        target_len,
                        target_start: b_beg,
                        target_end: b_end,
                        matches,
                        block_len,
                        mapping_quality: 60, // Default
                        cigar: String::new(), // Not extracting CIGAR for now
                        tags: Vec::new(),
                        mismatches: diffs - matches, // Rough approximation
                        gap_opens: 0, // Not tracked in .1aln
                        gap_len: 0,
                    };

                    return Ok(Some(alignment));
                }
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
