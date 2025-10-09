/// ONElib bindings for reading/writing .1aln files
///
/// This module provides safe Rust wrappers around the ONElib C library
/// for working with .1aln alignment files, using the onecode crate.

use std::path::Path;
use std::sync::Mutex;
use anyhow::{Result, bail, Context};
use onecode::{OneFile, OneSchema};

use crate::alignment::Alignment;

/// Global lock to serialize schema creation (prevents remaining race conditions in ONElib)
///
/// Despite mkstemp() fixes in onecode-rs, there are still race conditions when
/// multiple threads create schemas concurrently. This serializes creation as a workaround.
///
/// See: https://github.com/pangenome/onecode-rs/blob/main/SCHEMA_CLEANUP_BUG.md
static SCHEMA_CREATION_LOCK: Mutex<()> = Mutex::new(());

/// Create the schema for .1aln files (shared between reader and writer)
fn create_aln_schema() -> Result<OneSchema> {
    // Serialize schema creation to avoid race conditions
    let _guard = SCHEMA_CREATION_LOCK.lock().unwrap();

    let schema_text = r#"
P 3 aln
D t 1 3 INT
O g 0
G S 0
O S 1 6 STRING
D G 1 3 INT
D C 1 3 INT
D M 1 8 INT_LIST
O a 0
G A 0
D p 2 3 INT 3 INT
O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
D L 2 3 INT 3 INT
D R 0
D D 1 3 INT
D T 1 8 INT_LIST
D X 1 8 INT_LIST
D Q 1 3 INT
D E 1 3 INT
D Z 1 6 STRING
"#;

    OneSchema::from_text(schema_text)
        .context("Failed to create .1aln schema")
}

/// Alignment record with numeric IDs (compatible with C FFI API)
///
/// This struct maintains API compatibility with the old C FFI implementation
/// while using the new onecode-rs backend.
#[derive(Debug, Clone)]
pub struct AlnRecord {
    pub query_id: i64,
    pub query_start: i64,
    pub query_end: i64,
    pub target_id: i64,
    pub target_start: i64,
    pub target_end: i64,
    pub reverse: i32,
    pub diffs: i32,
    pub query_len: i64,
    pub target_len: i64,
}

/// Reader for .1aln files
pub struct AlnReader {
    file: OneFile,
    num_alignments: i64,
}

impl AlnReader {
    /// Open a .1aln file for reading
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_str = path.as_ref().to_str()
            .context("Invalid path")?;

        // Create schema for .1aln files
        let schema = create_aln_schema()?;

        let mut file = OneFile::open_read(path_str, Some(&schema), Some("aln"), 1)
            .context(format!("Failed to open .1aln file: {}", path_str))?;

        // Count alignments by scanning for 'A' records
        let num_alignments = Self::count_alignments(&mut file)?;

        Ok(AlnReader { file, num_alignments })
    }

    /// Count the total number of alignments in the file
    fn count_alignments(file: &mut OneFile) -> Result<i64> {
        let mut count = 0i64;

        loop {
            let line_type = file.read_line();
            if line_type == '\0' {
                break; // EOF
            }
            if line_type == 'A' {
                count += 1;
            }
        }

        // Reset to beginning
        file.goto('A', 0)?;

        Ok(count)
    }

    /// Get total number of alignments in the file
    pub fn num_alignments(&self) -> i64 {
        self.num_alignments
    }

    /// Read next alignment from the file
    /// Returns None when EOF is reached
    pub fn read_alignment(&mut self) -> Result<Option<Alignment>> {
        // Skip to next 'A' record (alignment)
        loop {
            let line_type = self.file.line_type();

            // If we're not at an 'A' record, read next line
            if line_type != 'A' {
                let next = self.file.read_line();
                if next == '\0' {
                    return Ok(None); // EOF
                }
                continue;
            }

            // We're at an 'A' record
            // Schema: O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
            // Fields: aread, abpos, aepos, bread, bbpos, bepos
            let a_id = self.file.int(0);
            let a_beg = self.file.int(1) as u64;
            let a_end = self.file.int(2) as u64;
            let b_id = self.file.int(3);
            let b_beg = self.file.int(4) as u64;
            let b_end = self.file.int(5) as u64;

            // Read associated data lines until we hit 'T' (trace) or next 'A'
            let mut matches = 0u64;
            let mut diffs = 0u64;
            let mut is_reverse = false;
            let mut query_len = 0u64;
            let mut target_len = 0u64;

            loop {
                let next_type = self.file.read_line();

                if next_type == '\0' {
                    break; // EOF
                }

                match next_type {
                    'T' => {
                        // Trace record - marks end of this alignment's data
                        // Read next line to position for next alignment
                        self.file.read_line();
                        break;
                    }
                    'M' => matches = self.file.int(0) as u64,
                    'D' => diffs = self.file.int(0) as u64,
                    'R' => is_reverse = true,
                    'L' => {
                        query_len = self.file.int(0) as u64;
                        target_len = self.file.int(1) as u64;
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

    /// Read all alignments from the file
    pub fn read_all(&mut self) -> Result<Vec<Alignment>> {
        let mut alignments = Vec::new();

        while let Some(alignment) = self.read_alignment()? {
            alignments.push(alignment);
        }

        Ok(alignments)
    }

    /// Get sequence name by ID (requires embedded GDB in .1aln file)
    ///
    /// Note: The `which_db` parameter is kept for API compatibility with the old
    /// C FFI implementation, but is ignored since onecode-rs automatically
    /// handles both query and target databases from the embedded GDB.
    pub fn get_seq_name(&mut self, seq_id: i64, _which_db: i32) -> Result<String> {
        self.file.get_sequence_name(seq_id)
            .ok_or_else(|| anyhow::anyhow!("Sequence ID {} not found", seq_id))
    }

    /// Get all sequence names at once (efficient for bulk lookups)
    pub fn get_all_seq_names(&mut self) -> std::collections::HashMap<i64, String> {
        self.file.get_all_sequence_names()
    }

    /// Read next alignment record (C FFI-compatible API)
    ///
    /// This method provides API compatibility with the old C FFI implementation.
    /// Returns `Ok(Some(record))` if a record was read, `Ok(None)` at EOF.
    pub fn read_record(&mut self) -> Result<Option<AlnRecord>> {
        // Read using the Alignment-based API
        let alignment = match self.read_alignment()? {
            Some(aln) => aln,
            None => return Ok(None),
        };

        // Convert Alignment to AlnRecord
        let query_id: i64 = alignment.query_name.parse()
            .unwrap_or(-1); // Default to -1 if not numeric
        let target_id: i64 = alignment.target_name.parse()
            .unwrap_or(-1);

        let record = AlnRecord {
            query_id,
            query_start: alignment.query_start as i64,
            query_end: alignment.query_end as i64,
            target_id,
            target_start: alignment.target_start as i64,
            target_end: alignment.target_end as i64,
            reverse: if alignment.strand == '-' { 1 } else { 0 },
            diffs: (alignment.mismatches + alignment.gap_len) as i32,
            query_len: alignment.query_len as i64,
            target_len: alignment.target_len as i64,
        };

        Ok(Some(record))
    }
}

/// Writer for .1aln files
pub struct AlnWriter {
    file: OneFile,
}

impl AlnWriter {
    /// Create a new .1aln file for writing
    pub fn create<P: AsRef<Path>>(path: P, binary: bool) -> Result<Self> {
        let path_str = path.as_ref().to_str()
            .context("Invalid path")?;

        // Create schema for .1aln files
        let schema = create_aln_schema()?;

        let file = OneFile::open_write_new(path_str, &schema, "aln", binary, 1)
            .context(format!("Failed to create .1aln file: {}", path_str))?;

        Ok(AlnWriter { file })
    }

    /// Write an alignment to the file
    pub fn write_alignment(&mut self, aln: &Alignment) -> Result<()> {
        // Parse sequence IDs from names (they should be numeric IDs for now)
        let a_id: i64 = aln.query_name.parse()
            .with_context(|| format!("Query name '{}' is not a numeric ID", aln.query_name))?;
        let b_id: i64 = aln.target_name.parse()
            .with_context(|| format!("Target name '{}' is not a numeric ID", aln.target_name))?;

        // Write 'A' line: alignment coordinates
        // Schema: O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT
        // Fields: aread, abpos, aepos, bread, bbpos, bepos
        self.file.set_int(0, a_id);
        self.file.set_int(1, aln.query_start as i64);
        self.file.set_int(2, aln.query_end as i64);
        self.file.set_int(3, b_id);
        self.file.set_int(4, aln.target_start as i64);
        self.file.set_int(5, aln.target_end as i64);
        self.file.write_line('A', 0, None);

        // Write 'R' line if reverse complement
        if aln.strand == '-' {
            self.file.write_line('R', 0, None);
        }

        // Write 'L' line: sequence lengths
        self.file.set_int(0, aln.query_len as i64);
        self.file.set_int(1, aln.target_len as i64);
        self.file.write_line('L', 0, None);

        // Write 'M' line: matches
        self.file.set_int(0, aln.matches as i64);
        self.file.write_line('M', 0, None);

        // Write 'D' line: differences (mismatches + gaps)
        let diffs = aln.mismatches + aln.gap_len;
        self.file.set_int(0, diffs as i64);
        self.file.write_line('D', 0, None);

        // Write 'Q' line: quality
        self.file.set_int(0, aln.mapping_quality as i64);
        self.file.write_line('Q', 0, None);

        // For now, we'll write empty trace lines
        // TODO: Generate actual trace points if needed
        self.file.write_line('T', 0, None);
        self.file.write_line('X', 0, None);

        Ok(())
    }

    /// Add provenance information to the file header
    pub fn add_provenance(&mut self, prog: &str, version: &str, command: &str) -> Result<()> {
        self.file.add_provenance(prog, version, command)?;
        Ok(())
    }

    /// Add reference file information to the header
    pub fn add_reference(&mut self, filename: &str, count: i64) -> Result<()> {
        self.file.add_reference(filename, count)?;
        Ok(())
    }

    /// Finalize the file (called automatically on drop)
    pub fn finalize(self) {
        self.file.close();
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
        assert!(!alignments.is_empty());
    }

    #[test]
    fn test_write_simple_alignment() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path();

        let mut writer = AlnWriter::create(temp_path, true).unwrap();

        let aln = Alignment {
            query_name: "0".to_string(),
            query_len: 10000,
            query_start: 1000,
            query_end: 1500,
            strand: '+',
            target_name: "1".to_string(),
            target_len: 20000,
            target_start: 2000,
            target_end: 2500,
            matches: 475,
            block_len: 500,
            mapping_quality: 60,
            cigar: String::new(),
            tags: Vec::new(),
            mismatches: 25,
            gap_opens: 0,
            gap_len: 0,
        };

        writer.write_alignment(&aln).unwrap();
        writer.finalize();
    }

    #[test]
    fn test_roundtrip() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path().to_path_buf();

        // Create a test alignment
        let original_aln = Alignment {
            query_name: "0".to_string(),
            query_len: 10000,
            query_start: 1000,
            query_end: 1500,
            strand: '+',
            target_name: "1".to_string(),
            target_len: 20000,
            target_start: 2000,
            target_end: 2500,
            matches: 475,
            block_len: 500,
            mapping_quality: 60,
            cigar: String::new(),
            tags: Vec::new(),
            mismatches: 25,
            gap_opens: 0,
            gap_len: 0,
        };

        // Write it
        {
            let mut writer = AlnWriter::create(&temp_path, true).unwrap();
            writer.write_alignment(&original_aln).unwrap();
            writer.finalize();
        }

        // Read it back
        let mut reader = AlnReader::open(&temp_path).unwrap();
        let alignments = reader.read_all().unwrap();

        // Verify
        assert_eq!(alignments.len(), 1, "Should read exactly one alignment");
        let read_aln = &alignments[0];

        assert_eq!(read_aln.query_name, original_aln.query_name);
        assert_eq!(read_aln.query_len, original_aln.query_len);
        assert_eq!(read_aln.query_start, original_aln.query_start);
        assert_eq!(read_aln.query_end, original_aln.query_end);
        assert_eq!(read_aln.strand, original_aln.strand);
        assert_eq!(read_aln.target_name, original_aln.target_name);
        assert_eq!(read_aln.target_len, original_aln.target_len);
        assert_eq!(read_aln.target_start, original_aln.target_start);
        assert_eq!(read_aln.target_end, original_aln.target_end);
        assert_eq!(read_aln.matches, original_aln.matches);
        assert_eq!(read_aln.mapping_quality, original_aln.mapping_quality);
    }
}
