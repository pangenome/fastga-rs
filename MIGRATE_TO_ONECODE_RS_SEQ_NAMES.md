# Migration Guide: Use onecode-rs for Sequence Name Extraction

## Overview

`onecode-rs` now supports sequence name extraction from embedded GDB in .1aln files. This means `fastga-rs` should consolidate to use `onecode-rs` for **all** ONElib operations, removing the custom FFI bindings in `src/aln_reader.rs`.

## Current State

`fastga-rs` has **two different** .1aln implementations:

1. **src/onelib.rs** - Uses `onecode-rs::OneFile` for reading/writing ✅
2. **src/aln_reader.rs** - Uses custom C FFI for sequence names ❌
   ```rust
   extern "C" {
       fn aln_get_seq_name(handle: *mut c_void, seq_id: i64, which_db: c_int) -> *const c_char;
   }
   ```

## What Changed in onecode-rs

New methods on `OneFile` (v0.1.0+):

```rust
impl OneFile {
    /// Get a single sequence name by ID
    pub fn get_sequence_name(&mut self, seq_id: i64) -> Option<String>

    /// Get all sequence names at once (efficient for bulk lookups)
    pub fn get_all_sequence_names(&mut self) -> HashMap<i64, String>
}
```

### Usage Example

```rust
use onecode::OneFile;

let mut file = OneFile::open_read("alignments.1aln", None, None, 1)?;

// Option 1: Get all names upfront (efficient)
let seq_names = file.get_all_sequence_names();

// Read alignments
loop {
    let line_type = file.read_line();
    if line_type == '\0' { break; }

    if line_type == 'A' {
        let query_id = file.int(0);
        let target_id = file.int(3);

        let query_name = seq_names.get(&query_id);
        let target_name = seq_names.get(&target_id);
    }
}

// Option 2: Get individual names on-demand
let name = file.get_sequence_name(5)?;
```

## Migration Steps for fastga-rs

### 1. Update Cargo.toml

Ensure you're using the latest onecode-rs:

```toml
[dependencies]
onecode = { git = "https://github.com/pangenome/onecode-rs", branch = "main" }
```

### 2. Remove Custom FFI (src/aln_reader.rs)

**Before:**
```rust
extern "C" {
    fn aln_get_seq_name(handle: *mut c_void, seq_id: i64, which_db: c_int) -> *const c_char;
}

impl AlnReader {
    pub fn get_seq_name(&self, seq_id: i64, which_db: i32) -> Option<String> {
        unsafe {
            let name_ptr = aln_get_seq_name(self.handle, seq_id, which_db as c_int);
            if name_ptr.is_null() {
                None
            } else {
                Some(CStr::from_ptr(name_ptr).to_string_lossy().into_owned())
            }
        }
    }
}
```

**After:**
```rust
// Remove the extern "C" block entirely

impl AlnReader {
    pub fn get_seq_name(&mut self, seq_id: i64) -> Option<String> {
        // AlnReader wraps OneFile internally
        self.file.get_sequence_name(seq_id)
    }

    pub fn get_all_seq_names(&mut self) -> HashMap<i64, String> {
        self.file.get_all_sequence_names()
    }
}
```

### 3. Update Internal Structure

If `AlnReader` doesn't already expose the underlying `OneFile`, add it:

```rust
pub struct AlnReader {
    file: OneFile,  // Make this accessible
    // ... other fields
}

impl AlnReader {
    // Provide direct access if needed
    pub fn file_mut(&mut self) -> &mut OneFile {
        &mut self.file
    }
}
```

### 4. Update Tests

Replace FFI-based tests with onecode-rs-based tests:

```rust
#[test]
fn test_sequence_names() {
    let mut reader = AlnReader::open("test.1aln").unwrap();

    // Use the new consolidated API
    let names = reader.get_all_seq_names();
    assert!(!names.is_empty());

    let name0 = reader.get_seq_name(0).unwrap();
    assert!(name0.contains("expected_string"));
}
```

## Benefits of This Migration

✅ **Single source of truth**: All ONElib C bindings in one place (onecode-rs)
✅ **Less code duplication**: No need for multiple FFI implementations
✅ **Easier maintenance**: Update onecode-rs once, everyone benefits
✅ **Better type safety**: onecode-rs provides the safe Rust API layer
✅ **Consistent API**: Same patterns across all tools

## Impact on sweepga and Other Consumers

After this migration, `sweepga` (and other tools) can choose:

1. **Use onecode-rs directly** for reading .1aln with names
2. **Use fastga-rs's AlnReader** (which now uses onecode-rs internally)

Both approaches work correctly and use the same underlying implementation.

## Testing

After migration, verify:

1. All existing tests pass
2. Sequence name extraction works correctly
3. No performance regression (the new API should be as fast or faster)

```bash
cd ~/fastga-rs
cargo update  # Update onecode-rs dependency
cargo test
```

## Questions?

- See onecode-rs examples: `~/onecode-rs/tests/sequence_names_test.rs`
- Check onecode-rs documentation: `~/onecode-rs/README.md`
- The implementation uses standard ONElib operations (oneGoto, oneReadLine) to scan S lines

## Implementation Notes

The onecode-rs implementation:
- Scans for `S` (scaffold/sequence) lines which contain embedded GDB data
- Sequence IDs are 0-indexed and correspond to S line order
- Uses `oneGoto()` to navigate efficiently through binary files
- Preserves current file position after lookups
- Returns `HashMap<i64, String>` for bulk operations (faster than repeated lookups)
