//! Native FFI bindings to modified FastGA C code for true streaming.
//!
//! This module provides direct integration with FastGA's alignment generation,
//! intercepting alignments as they are created rather than parsing output files.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void, c_double};
use std::path::Path;
use std::sync::Mutex;
use crate::alignment::Alignment;
use crate::error::{Result, FastGAError};

/// Represents a FastGA overlap structure
#[repr(C)]
pub struct FastGAOverlap {
    pub path: FastGAPath,
    pub flags: u32,
    pub aread: c_int,
    pub bread: c_int,
}

/// Represents a FastGA path/alignment
#[repr(C)]
pub struct FastGAPath {
    pub trace: *mut c_void,
    pub tlen: c_int,
    pub diffs: c_int,
    pub abpos: c_int,
    pub bbpos: c_int,
    pub aepos: c_int,
    pub bepos: c_int,
}

/// C callback function type
type FastGACallback = unsafe extern "C" fn(
    user_data: *mut c_void,
    ovl: *const FastGAOverlap,
    aln: *mut c_void,  // Alignment structure
    query_name: *const c_char,
    target_name: *const c_char,
) -> c_int;

// External C functions from our modified FastGA
extern "C" {
    fn fastga_align_with_callback(
        genome1_path: *const c_char,
        genome2_path: *const c_char,
        callback: FastGACallback,
        user_data: *mut c_void,
        num_threads: c_int,
        min_length: c_int,
        min_identity: c_double,
    ) -> c_int;
}

/// Container for the Rust callback
struct CallbackContainer<F>
where
    F: FnMut(Alignment) -> bool
{
    callback: F,
    processed: usize,
    kept: usize,
}

/// Convert FastGA overlap to Rust alignment
unsafe fn overlap_to_alignment(
    ovl: *const FastGAOverlap,
    query_name: *const c_char,
    target_name: *const c_char,
) -> Result<Alignment> {
    let ovl = &*ovl;
    let path = &ovl.path;

    let query = CStr::from_ptr(query_name).to_string_lossy().into_owned();
    let target = CStr::from_ptr(target_name).to_string_lossy().into_owned();

    // Calculate alignment metrics
    let query_start = path.abpos as usize;
    let query_end = path.aepos as usize;
    let target_start = path.bbpos as usize;
    let target_end = path.bepos as usize;

    let query_len = query_end - query_start;
    let target_len = target_end - target_start;

    // Simplified - in reality we'd decode the trace to get exact matches/mismatches
    let matches = ((query_len.min(target_len) as f64) * 0.9) as usize;
    let mismatches = path.diffs as usize;

    let strand = if (ovl.flags & 0x1) != 0 { '-' } else { '+' };

    Ok(Alignment {
        query_name: query,
        query_len: query_end,  // This should be total sequence length
        query_start,
        query_end,
        strand,
        target_name: target,
        target_len: target_end,  // This should be total sequence length
        target_start,
        target_end,
        matches,
        block_len: query_len,
        mapping_quality: 255,
        cigar: generate_cigar_from_trace(path),
        mismatches,
        gap_opens: 0,
        gap_len: 0,
        tags: Vec::new(),
    })
}

/// Generate CIGAR string from trace points
fn generate_cigar_from_trace(path: &FastGAPath) -> String {
    // Simplified CIGAR generation
    // In production, we'd decode the actual trace points
    if path.diffs == 0 {
        format!("{}=", path.aepos - path.abpos)
    } else {
        // Mix of matches and mismatches
        let total = path.aepos - path.abpos;
        let matches = total - path.diffs;
        format!("{}={}X", matches, path.diffs)
    }
}

/// C callback wrapper that forwards to Rust closure
unsafe extern "C" fn callback_wrapper<F>(
    user_data: *mut c_void,
    ovl: *const FastGAOverlap,
    _aln: *mut c_void,
    query_name: *const c_char,
    target_name: *const c_char,
) -> c_int
where
    F: FnMut(Alignment) -> bool
{
    let container = &mut *(user_data as *mut Mutex<CallbackContainer<F>>);
    let mut container = container.lock().unwrap();

    container.processed += 1;

    match overlap_to_alignment(ovl, query_name, target_name) {
        Ok(alignment) => {
            if (container.callback)(alignment) {
                container.kept += 1;
                1  // Keep alignment
            } else {
                0  // Skip alignment
            }
        }
        Err(e) => {
            eprintln!("Error converting alignment: {}", e);
            1  // Continue on error
        }
    }
}

/// Perform alignment with native FastGA integration
pub fn align_native<F>(
    genome1: &Path,
    genome2: &Path,
    num_threads: usize,
    min_length: usize,
    min_identity: f64,
    callback: F,
) -> Result<(usize, usize)>
where
    F: FnMut(Alignment) -> bool,
{
    let genome1_str = genome1.to_str()
        .ok_or_else(|| FastGAError::Other("Invalid genome1 path".to_string()))?;
    let genome2_str = genome2.to_str()
        .ok_or_else(|| FastGAError::Other("Invalid genome2 path".to_string()))?;

    let genome1_c = CString::new(genome1_str)?;
    let genome2_c = CString::new(genome2_str)?;

    let container = Mutex::new(CallbackContainer {
        callback,
        processed: 0,
        kept: 0,
    });

    unsafe {
        let result = fastga_align_with_callback(
            genome1_c.as_ptr(),
            genome2_c.as_ptr(),
            callback_wrapper::<F>,
            &container as *const _ as *mut c_void,
            num_threads as c_int,
            min_length as c_int,
            min_identity,
        );

        if result < 0 {
            return Err(FastGAError::FfiError("FastGA alignment failed".to_string()));
        }
    }

    let final_container = container.into_inner().unwrap();
    Ok((final_container.processed, final_container.kept))
}