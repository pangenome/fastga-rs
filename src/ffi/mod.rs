//! Foreign Function Interface (FFI) bindings to FastGA C code.
//!
//! This module provides low-level unsafe bindings to the FastGA C library.
//! These bindings are wrapped by the safe API in the parent module.

pub mod streaming;

use libc::{c_char, c_int, c_void, size_t};
use std::ffi::CString;
use std::ptr;

/// Opaque type representing a GDB (Genome Database)
#[repr(C)]
pub struct GDB {
    _private: [u8; 0],
}

/// Opaque type representing an alignment
#[repr(C)]
pub struct Alignment {
    _private: [u8; 0],
}

/// Opaque type representing a path/trace
#[repr(C)]
pub struct Path {
    pub trace: *mut c_void,
    pub tlen: c_int,
    pub diffs: c_int,
    pub abpos: c_int,
    pub bbpos: c_int,
    pub aepos: c_int,
    pub bepos: c_int,
}

/// Opaque type representing an overlap
#[repr(C)]
pub struct Overlap {
    pub path: Path,
    pub flags: u32,
    pub aread: c_int,
    pub bread: c_int,
}

// FFI function declarations would go here
// For now, these are placeholder declarations that would need to be
// properly linked to the compiled FastGA C code

extern "C" {
    // GDB functions
    pub fn Read_GDB(gdb: *mut GDB, fname: *const c_char) -> c_int;
    pub fn Close_GDB(gdb: *mut GDB);

    // Alignment functions (placeholder - actual signatures TBD)
    // pub fn Local_Alignment(...) -> *mut Alignment;
    // pub fn Compute_Trace(...) -> c_int;
}

/// Safe wrapper for GDB operations
pub struct SafeGDB {
    ptr: *mut GDB,
}

impl SafeGDB {
    /// Creates a new GDB from a file
    pub fn from_file(path: &str) -> Result<Self, String> {
        unsafe {
            let c_path = CString::new(path).map_err(|e| e.to_string())?;
            let gdb = libc::malloc(std::mem::size_of::<GDB>()) as *mut GDB;
            if gdb.is_null() {
                return Err("Failed to allocate GDB".to_string());
            }

            let result = Read_GDB(gdb, c_path.as_ptr());
            if result != 0 {
                libc::free(gdb as *mut c_void);
                return Err(format!("Failed to read GDB: error code {}", result));
            }

            Ok(SafeGDB { ptr: gdb })
        }
    }
}

impl Drop for SafeGDB {
    fn drop(&mut self) {
        unsafe {
            if !self.ptr.is_null() {
                Close_GDB(self.ptr);
                libc::free(self.ptr as *mut c_void);
            }
        }
    }
}

// Additional FFI bindings will be added as we develop the streaming API