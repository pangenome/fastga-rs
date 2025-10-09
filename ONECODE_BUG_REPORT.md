# Bug Report for ONEcode: Thread-Safety Issue in oneSchemaCreateFromText()

**Upstream Issue**: https://github.com/thegenemyers/ONEcode/issues/5 ✅ **FILED**
**Repository**: https://github.com/thegenemyers/ONEcode
**File**: `ONElib.c`, lines 394-421
**Severity**: High (causes FATAL ERROR in multi-threaded applications)

---

## Problem

`oneSchemaCreateFromText()` uses a single static temp file path based on PID, causing race conditions when multiple threads create schemas concurrently.

## Affected Code

```c
OneSchema *oneSchemaCreateFromText (const char *text)
{
  static char template[64] ;
  sprintf (template, "/tmp/OneTextSchema-%d.schema", getpid()) ;  // ← All threads use SAME path!

  errno = 0 ;
  FILE *f = fopen (template, "w") ;
  // ... write schema ...
  fclose (f) ;

  OneSchema *vs = oneSchemaCreateFromFile (template) ;

  errno = 0 ;
  unlink (template) ;  // ← Race condition here!
  if (errno) die ("failed to remove temporary file %s errno %d\n", template, errno) ;

  return vs ;
}
```

## Symptoms

When multiple threads call `oneSchemaCreateFromText()` in parallel:

```
FATAL ERROR: failed to remove temporary file /tmp/OneTextSchema-12345.schema errno 2
```

Or corrupted schema reads when threads overwrite each other's temp files.

## Reproduction

**Test code** (Rust, but issue is in C library):

```rust
#[test]
fn test_parallel_schema_creation_1() {
    let schema = OneSchema::from_text("P 3 aln\nO A 3 3 INT 3 INT 3 INT\n").unwrap();
    // ... use schema ...
}

#[test]
fn test_parallel_schema_creation_2() {
    let schema = OneSchema::from_text("P 3 aln\nO A 3 3 INT 3 INT 3 INT\n").unwrap();
    // ... use schema ...
}
```

**Run:**
```bash
cargo test  # ❌ FAILS: Both tests use /tmp/OneTextSchema-<pid>.schema
cargo test -- --test-threads=1  # ✅ Works with serial execution
```

## Root Cause

1. **Static buffer reuse**: `static char template[64]` is shared across all calls
2. **PID-only uniqueness**: `/tmp/OneTextSchema-<pid>.schema` is the same for all threads in a process
3. **Race condition in cleanup**: Multiple threads try to unlink the same file

## Proposed Fixes

### Option 1: Use `mkstemp()` for guaranteed unique temp files ✅ RECOMMENDED

```c
OneSchema *oneSchemaCreateFromText (const char *text)
{
  char template[] = "/tmp/OneTextSchema-XXXXXX";

  errno = 0;
  int fd = mkstemp(template);  // Creates unique file atomically
  if (fd < 0) die ("failed to create temporary file for schema");

  FILE *f = fdopen(fd, "w");
  if (!f) {
    close(fd);
    unlink(template);
    die ("failed to open temporary file %s for writing schema", template);
  }

  // ... write schema ...
  fclose (f);

  OneSchema *vs = oneSchemaCreateFromFile (template);

  errno = 0;
  unlink (template);
  if (errno) die ("failed to remove temporary file %s errno %d\n", template, errno);

  return vs;
}
```

**Benefits:**
- Thread-safe by design (atomic unique file creation)
- No race conditions
- POSIX standard, widely supported
- Minimal code changes

### Option 2: Add thread ID to filename

```c
#include <pthread.h>

OneSchema *oneSchemaCreateFromText (const char *text)
{
  char template[128];
  sprintf(template, "/tmp/OneTextSchema-%d-%lu.schema",
          getpid(), (unsigned long)pthread_self());
  // ... rest of implementation ...
}
```

**Drawbacks:**
- Requires pthread dependency
- Less clean than mkstemp()
- Still has potential (unlikely) collisions if thread IDs wrap

### Option 3: Use atomic counter with PID

```c
#include <stdatomic.h>

OneSchema *oneSchemaCreateFromText (const char *text)
{
  static atomic_int counter = 0;
  char template[128];
  sprintf(template, "/tmp/OneTextSchema-%d-%d.schema",
          getpid(), atomic_fetch_add(&counter, 1));
  // ... rest of implementation ...
}
```

**Drawbacks:**
- Requires C11 atomics
- More complex than mkstemp()

## Workaround (Currently Used)

In Rust bindings (onecode-rs/fastga-rs), we serialize schema creation:

```rust
static SCHEMA_CREATION_LOCK: Mutex<()> = Mutex::new(());

fn create_aln_schema() -> Result<OneSchema> {
    let _guard = SCHEMA_CREATION_LOCK.lock().unwrap();
    OneSchema::from_text(schema_text)  // Only one thread at a time
}
```

This works but prevents concurrent schema creation.

## Impact

- **fastga-rs**: Tests fail without workaround (serialization mutex)
- **onecode-rs**: Parallel test execution broken
- **Any multi-threaded app**: Risk of schema creation failures

## Environment

- ONElib version: Latest from https://github.com/thegenemyers/ONEcode
- Observed in: Rust bindings (onecode-rs), but issue is in C library
- Platform: Linux (all versions), macOS (likely)

## References

- Bug report: https://github.com/pangenome/onecode-rs/blob/main/SCHEMA_CLEANUP_BUG.md
- Workaround commit: https://github.com/pangenome/fastga-rs/commit/50fb1b1

## Recommendation

**Please implement Option 1 (`mkstemp()`)**. It's the cleanest, most thread-safe solution with minimal code changes and no new dependencies (POSIX standard).

Thank you for maintaining ONEcode! Happy to test any fixes.

---
**Filed by**: fastga-rs/onecode-rs maintainers
**Date**: 2025-10-09
