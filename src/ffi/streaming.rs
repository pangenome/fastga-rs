//! FFI bindings for streaming alignment processing.
//!
//! This module provides low-level unsafe bindings for streaming alignments
//! directly from FastGA without intermediate files.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void, c_double};
use std::ptr;
use std::sync::{Arc, Mutex};
use crate::alignment::Alignment;
use crate::error::{Result, FastGAError};

/// Raw alignment data passed from C callback
#[repr(C)]
pub struct RawAlignment {
    pub query_name: *const c_char,
    pub query_len: c_int,
    pub query_start: c_int,
    pub query_end: c_int,
    pub target_name: *const c_char,
    pub target_len: c_int,
    pub target_start: c_int,
    pub target_end: c_int,
    pub strand: c_int,
    pub cigar: *const c_char,
    pub matches: c_int,
    pub mismatches: c_int,
    pub gaps: c_int,
}

/// Type for the C callback function
pub type AlignmentCallbackC = unsafe extern "C" fn(
    user_data: *mut c_void,
    query_name: *const c_char,
    query_len: c_int,
    query_start: c_int,
    query_end: c_int,
    target_name: *const c_char,
    target_len: c_int,
    target_start: c_int,
    target_end: c_int,
    strand: c_int,
    cigar: *const c_char,
    matches: c_int,
    mismatches: c_int,
    gaps: c_int,
) -> c_int;

// External C functions
extern "C" {
    fn fastga_align_streaming(
        genome1_path: *const c_char,
        genome2_path: *const c_char,
        callback: AlignmentCallbackC,
        user_data: *mut c_void,
        num_threads: c_int,
        min_length: c_int,
        min_identity: c_double,
    ) -> c_int;
}

/// Context for streaming callbacks
pub struct StreamingContext<F>
where
    F: FnMut(Alignment) -> bool
{
    pub callback: F,
    pub alignments_processed: usize,
    pub alignments_kept: usize,
    pub error: Option<String>,
}

/// Convert raw C alignment data to Rust Alignment
unsafe fn raw_to_alignment(raw: &RawAlignment) -> Result<Alignment> {
    let query_name = CStr::from_ptr(raw.query_name)
        .to_string_lossy()
        .into_owned();

    let target_name = CStr::from_ptr(raw.target_name)
        .to_string_lossy()
        .into_owned();

    let cigar = CStr::from_ptr(raw.cigar)
        .to_string_lossy()
        .into_owned();

    let strand = if raw.strand == 0 { '+' } else { '-' };

    Ok(Alignment {
        query_name,
        query_len: raw.query_len as usize,
        query_start: raw.query_start as usize,
        query_end: raw.query_end as usize,
        strand,
        target_name,
        target_len: raw.target_len as usize,
        target_start: raw.target_start as usize,
        target_end: raw.target_end as usize,
        matches: raw.matches as usize,
        block_len: (raw.query_end - raw.query_start) as usize,
        mapping_quality: 255,  // Default high quality
        cigar,
        mismatches: raw.mismatches as usize,
        gap_opens: 0,  // Would need to parse from CIGAR
        gap_len: raw.gaps as usize,
        tags: Vec::new(),
    })
}

/// C callback that forwards to Rust closure
unsafe extern "C" fn alignment_callback_wrapper<F>(
    user_data: *mut c_void,
    query_name: *const c_char,
    query_len: c_int,
    query_start: c_int,
    query_end: c_int,
    target_name: *const c_char,
    target_len: c_int,
    target_start: c_int,
    target_end: c_int,
    strand: c_int,
    cigar: *const c_char,
    matches: c_int,
    mismatches: c_int,
    gaps: c_int,
) -> c_int
where
    F: FnMut(Alignment) -> bool
{
    let ctx = &mut *(user_data as *mut StreamingContext<F>);

    ctx.alignments_processed += 1;

    // Create raw alignment struct
    let raw = RawAlignment {
        query_name,
        query_len,
        query_start,
        query_end,
        target_name,
        target_len,
        target_start,
        target_end,
        strand,
        cigar,
        matches,
        mismatches,
        gaps,
    };

    // Convert to Rust alignment
    match raw_to_alignment(&raw) {
        Ok(alignment) => {
            // Call user callback
            let keep = (ctx.callback)(alignment);
            if keep {
                ctx.alignments_kept += 1;
                1  // Continue and keep alignment
            } else {
                0  // Skip this alignment
            }
        }
        Err(e) => {
            ctx.error = Some(e.to_string());
            -1  // Stop processing on error
        }
    }
}

/// Safe wrapper for streaming alignment
pub fn align_streaming<F>(
    genome1: &str,
    genome2: &str,
    num_threads: usize,
    min_length: usize,
    min_identity: f64,
    mut callback: F,
) -> Result<(usize, usize)>
where
    F: FnMut(Alignment) -> bool,
{
    let genome1_c = CString::new(genome1)
        .map_err(|e| FastGAError::Other(format!("Invalid genome1 path: {}", e)))?;
    let genome2_c = CString::new(genome2)
        .map_err(|e| FastGAError::Other(format!("Invalid genome2 path: {}", e)))?;

    let mut ctx = StreamingContext {
        callback,
        alignments_processed: 0,
        alignments_kept: 0,
        error: None,
    };

    unsafe {
        let result = fastga_align_streaming(
            genome1_c.as_ptr(),
            genome2_c.as_ptr(),
            alignment_callback_wrapper::<F>,
            &mut ctx as *mut _ as *mut c_void,
            num_threads as c_int,
            min_length as c_int,
            min_identity as c_double,
        );

        if result < 0 {
            if let Some(error) = ctx.error {
                return Err(FastGAError::FfiError(error));
            } else {
                return Err(FastGAError::FfiError("FastGA streaming failed".to_string()));
            }
        }
    }

    Ok((ctx.alignments_processed, ctx.alignments_kept))
}