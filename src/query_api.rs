/// Query-specific API for FastGA alignments
/// Ensures complete alignment results for each query sequence

use std::ffi::{CString, CStr};
use std::path::Path;
use std::os::raw::{c_char, c_int, c_void};
use anyhow::Result;

use crate::{Alignment, AlignmentSet};
use crate::error::FastGAError;

/// Represents a single query sequence to be aligned
pub struct Query {
    pub name: String,
    pub sequence: Vec<u8>,
}

/// Results from aligning a single query
pub struct QueryResult {
    pub query_name: String,
    pub alignments: Vec<Alignment>,
}

/// Query-specific aligner that guarantees complete results per query
pub struct QueryAligner {
    target_gdb: String,  // Path to target GDB
    config: crate::Config,
}

impl QueryAligner {
    /// Create a new query-specific aligner for a target database
    pub fn new(target_db: impl AsRef<Path>, config: crate::Config) -> Result<Self> {
        // Convert target to GDB if needed
        let target_gdb = crate::ffi::prepare_database(target_db.as_ref())?;

        Ok(Self {
            target_gdb,
            config,
        })
    }

    /// Align a single query and get all its alignments
    /// This guarantees you get ALL alignments for this query before returning
    pub fn align_query(&self, query: &Query) -> Result<QueryResult> {
        // Create temporary GDB for single query
        let temp_dir = tempfile::tempdir()?;
        let query_fasta = temp_dir.path().join("query.fa");

        // Write query to FASTA
        std::fs::write(&query_fasta, format!(
            ">{}\n{}\n",
            query.name,
            std::str::from_utf8(&query.sequence)?
        ))?;

        // Convert to GDB
        let query_gdb = crate::ffi::prepare_database(&query_fasta)?;

        // Run alignment for just this query
        let alignments = crate::ffi::run_fastga_alignment(
            &query_gdb,
            &self.target_gdb,
            &self.config
        )?;

        // Filter to only this query's alignments (should be all of them)
        let query_alignments: Vec<Alignment> = alignments.alignments
            .into_iter()
            .filter(|a| a.query_name == query.name)
            .collect();

        Ok(QueryResult {
            query_name: query.name.clone(),
            alignments: query_alignments,
        })
    }

    /// Align multiple queries, processing each completely before moving to the next
    /// The callback is called once per query with all its alignments
    pub fn align_queries<F>(&self, queries: Vec<Query>, mut callback: F) -> Result<()>
    where
        F: FnMut(QueryResult) -> bool,
    {
        for query in queries {
            let result = self.align_query(&query)?;

            // Callback gets complete results for this query
            if !callback(result) {
                break;  // Stop if callback returns false
            }
        }

        Ok(())
    }

    /// Stream alignments but with query boundaries marked
    /// Calls start_query/end_query to delimit query boundaries
    pub fn align_queries_with_boundaries<F, G, H>(
        &self,
        queries: Vec<Query>,
        mut on_query_start: F,
        mut on_alignment: G,
        mut on_query_end: H,
    ) -> Result<()>
    where
        F: FnMut(&str),           // Called when starting a query
        G: FnMut(&Alignment),     // Called for each alignment
        H: FnMut(&str, usize),    // Called when query is complete with count
    {
        for query in queries {
            on_query_start(&query.name);

            let result = self.align_query(&query)?;
            let count = result.alignments.len();

            for alignment in &result.alignments {
                on_alignment(alignment);
            }

            on_query_end(&query.name, count);
        }

        Ok(())
    }
}

/// Alternative: Patch FastGA to use a modified version that processes queries sequentially
/// This would be more efficient than creating temporary GDBs
pub mod ffi_patch {
    use super::*;

    // Link to our patched C API
    extern "C" {
        // Our custom function that processes one query at a time
        fn align_single_query(
            query_gdb: *const c_void,
            query_idx: c_int,
            target_gdb: *const c_void,
            work: *mut c_void,
            spec: *mut c_void,
            callback: extern "C" fn(*const c_void, *mut c_void) -> c_int,
            user_data: *mut c_void,
        ) -> *mut c_void;

        fn free_query_alignments(set: *mut c_void);
    }

    /// Low-level FFI wrapper for single-query alignment
    pub unsafe fn align_query_ffi(
        query_gdb: *const c_void,
        query_idx: i32,
        target_gdb: *const c_void,
    ) -> Vec<Alignment> {
        // This would call our C patch
        // For now, this is a placeholder showing how it would work
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query_alignment() {
        // Test that we get complete results for each query
        let query = Query {
            name: "test_query".to_string(),
            sequence: b"ACGTACGTACGT".to_vec(),
        };

        // Would test against a target database
        // Verify all alignments for this query are returned
    }

    #[test]
    fn test_query_boundaries() {
        // Test that query boundaries are properly maintained
        // Each query should be fully processed before the next
    }
}