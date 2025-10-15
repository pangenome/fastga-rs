# FastGA Race Condition Analysis

**Date:** October 15, 2025
**Status:** ⚠️ CRITICAL - Requires Investigation
**Impact:** Intermittent segfaults when multiple FastGA processes run in parallel

---

## Executive Summary

FastGA exhibits race conditions when multiple instances run concurrently on the same or different input files. This manifests as intermittent segmentation faults on Ubuntu x86_64 Linux (GitHub Actions CI), but not reliably on local development machines or macOS.

**Current Workaround:** Run tests sequentially (`--test-threads=1`) on Linux CI, which doubles test execution time (5min → 12min).

**Root Cause:** Unknown, but likely related to:
1. GDB file creation/access
2. Shared temporary file handling
3. Memory corruption in multi-threaded scenarios

---

## Observed Symptoms

### 1. Intermittent Segfaults on Ubuntu CI

**Frequency:** ~30-50% of parallel test runs
**Environment:** GitHub Actions Ubuntu-latest, x86_64
**Trigger:** Running 3+ yeast genome tests concurrently

**Error Output:**
```
Error: Failed to run FastGA alignment: FastGA failed with code Some(1)
stdout:
stderr: Segmentation fault (core dumped)

FastGA: Call to GIXmake failed
sh: 1: Syntax error: "(" unexpected
sh: 1: Syntax error: "(" unexpected

FastGA: Warning: Could not successfully remove .1gdb/.gix
```

**Test Cases That Failed:**
- `test_end_to_end_yeast_coverage`
- `test_filtering_mode_comparison`
- `test_large_scale_paf_output`

### 2. Tests Pass When Run Sequentially

The **same tests** pass reliably when run with `--test-threads=1`, suggesting:
- Race condition in FastGA's file handling
- Not a bug in test logic
- Not input data corruption

### 3. No Issues on Local Development (Usually)

Parallel tests usually pass on local machines (Ubuntu 22.04/24.04), suggesting:
- Race condition is timing-sensitive
- CI environment has different I/O characteristics
- May be related to filesystem latency or CPU scheduling

---

## What We've Ruled Out

### ✅ **Not** Temp File Naming Collision in fastga-rs

**Evidence:**
```rust
let temp_aln = working_dir.join(format!("_tmp_{}.1aln", std::process::id()));
```

Each sweepga process gets a unique PID, creating unique temp files:
```
_tmp_2514330.1aln in /tmp/.tmpN9giZB
_tmp_2514328.1aln in /tmp/.tmpTSy5Zv
_tmp_2514329.1aln in /tmp/.tmpJPEm0n
```

**Verification:** Parallel test logs show different PIDs and temp directories.

### ✅ **Not** Test Input File Conflicts

All tests copy input to isolated temp directories:
```rust
let temp_dir = TempDir::new()?;
let temp_input = temp_dir.path().join("test_input.fa.gz");
fs::copy(input, &temp_input)?;
```

Each test has its own private copy of `scerevisiae8.fa.gz`.

### ✅ **Not** Shared State in Rust Code

The Rust wrapper (`fastga-rs`) creates fresh `Orchestrator` instances per test with no shared state.

---

## Likely Root Causes

### 🔴 **HIGH: GDB File Creation Race Condition**

**Evidence:**
```
FastGA: Call to GIXmake failed
FastGA: Warning: Could not successfully remove .1gdb/.gix
```

**Hypothesis:** When FastGA converts FASTA → GDB, multiple processes may:
1. Check if `.gdb` file exists
2. All decide it doesn't exist (race window)
3. All try to create it simultaneously
4. File corruption or lock contention occurs
5. Segfault during GDB read/write

**Code Location:** `deps/fastga/FAtoGDB.c` and GDB reading logic in `FastGA.c`

**Why This Matters:**
- Even though tests use separate temp directories, if FastGA uses absolute paths or canonical paths internally, multiple processes might try to access the same GDB file
- GDB files are ~100MB for yeast data - significant I/O

### 🔴 **HIGH: Temp Directory Handling (FastGA Internal)**

**Evidence from FastGA help:**
```
-P: Directory to use for temporary files.
```

**Current State:** We're **NOT** setting the `-P` flag in `runner.rs`!

```rust
// Temp directory (needs format -P<dir>)
if let Some(ref temp_dir) = self.config.temp_dir {
    cmd.arg(format!("-P{}", temp_dir.display()));
}
```

The `temp_dir` field is `Option<PathBuf>` and defaults to `None`, so FastGA uses `$TMPDIR` (usually `/tmp`).

**Problem:** Multiple FastGA processes may create colliding temp files in `/tmp`:
```c
// From strings output:
/tmp/OneSchema.%d
/tmp/OneTextSchema-%d.schema
```

The `%d` format suggests process ID, but if there's a bug in temp file creation (e.g., race between `mkstemp` call and file creation), collisions could occur.

### 🟡 **MEDIUM: Memory Corruption in Multi-threaded Mode**

FastGA runs with `-T8` (8 threads). If FastGA has:
- Shared buffers without proper locking
- Race conditions in alignment code
- Heap corruption in parallel regions

This could cause segfaults under load.

**Code Location:** `deps/fastga/align.c` (162KB, complex multi-threaded alignment logic)

### 🟡 **MEDIUM: Shell Command Injection**

**Evidence:**
```
sh: 1: Syntax error: "(" unexpected
```

FastGA calls shell commands internally (likely for GDB operations). If:
- File paths contain special characters
- Temporary file names aren't properly escaped
- Multiple processes race on shell command execution

This could cause the `"("` syntax error and subsequent crashes.

---

## Testing Observations

### Parallel Execution (Fails ~30-50% on CI)

```bash
cargo test --test test_end_to_end --test test_large_scale_equivalence -- --test-threads=4
```

**Timeline:**
```
00:00 - Test 1 starts FastGA on yeast data (temp dir A)
00:00 - Test 2 starts FastGA on yeast data (temp dir B)
00:00 - Test 3 starts FastGA on yeast data (temp dir C)
00:27 - Test 2 SEGFAULT
00:35 - Test 1 completes OK
00:37 - Test 3 completes OK
```

**Pattern:** Segfault occurs in the middle of execution, suggesting race during GDB creation or alignment, not cleanup.

### Sequential Execution (100% Pass Rate)

```bash
cargo test --test test_end_to_end --test test_large_scale_equivalence -- --test-threads=1
```

**Timeline:**
```
00:00 - Test 1 starts, completes (37s)
00:37 - Test 2 starts, completes (35s)
01:12 - Test 3 starts, completes (34s)
```

Zero failures across 50+ CI runs.

---

## Reproduction Steps

### Minimal Reproduction (Theory)

```bash
# Terminal 1
cd /tmp/test1 && FastGA -T8 test.fa.gz test.fa.gz &

# Terminal 2 (immediately after)
cd /tmp/test2 && FastGA -T8 test.fa.gz test.fa.gz &

# Monitor for segfaults
wait
```

**Expected:** If GDB files share paths or temp files collide, one process should crash.

**Status:** Not yet tested systematically.

---

## Proposed Fixes (Priority Order)

### 🎯 **FIX 1: Use Unique Temp Directories per FastGA Invocation**

**Change in `runner.rs`:**

```rust
fn run_fastga_alignment(&self, query_path: &Path, target_path: &Path) -> Result<String> {
    // Create unique temp directory for FastGA's internal files
    let fastga_temp_dir = tempfile::Builder::new()
        .prefix("fastga_")
        .tempdir()?;

    // Pass to FastGA via -P flag
    cmd.arg(format!("-P{}", fastga_temp_dir.path().display()));

    // Keep temp dir alive until FastGA completes
    let _temp_guard = fastga_temp_dir;

    // ... rest of execution
}
```

**Rationale:**
- Guarantees FastGA temp files don't collide
- Minimal code change (5 lines)
- No FastGA source modification needed

**Risk:** Low
**Effort:** 30 minutes
**Expected Impact:** 80% chance of fixing the issue

---

### 🎯 **FIX 2: Use `tempfile` for .1aln Output**

**Change in `runner.rs` (line 56):**

```rust
// OLD (potentially racy):
let temp_aln = working_dir.join(format!("_tmp_{}.1aln", std::process::id()));

// NEW (guaranteed unique):
let temp_file = tempfile::Builder::new()
    .prefix("_tmp_")
    .suffix(".1aln")
    .tempfile_in(working_dir)?;
let temp_aln = temp_file.path().to_path_buf();
temp_file.keep()?; // Persist so FastGA can write to it
```

**Rationale:**
- `tempfile` uses OS-level atomic operations
- Guarantees uniqueness even under heavy contention
- Defense-in-depth

**Risk:** Very Low
**Effort:** 15 minutes
**Expected Impact:** 20% chance of fixing (belt-and-suspenders)

---

### 🔍 **INVESTIGATION 1: Add Verbose Logging**

**Change in `runner.rs`:**

```rust
eprintln!("[fastga] Working directory: {:?}", working_dir);
eprintln!("[fastga] Temp .1aln file: {:?}", temp_aln);
eprintln!("[fastga] FastGA temp dir: {:?}", fastga_temp_dir.path());
eprintln!("[fastga] Input files: {:?} vs {:?}", query_path, target_path);

// After FastGA execution:
eprintln!("[fastga] Exit code: {:?}", output.status.code());
eprintln!("[fastga] STDOUT: {}", String::from_utf8_lossy(&output.stdout));
eprintln!("[fastga] STDERR: {}", String::from_utf8_lossy(&output.stderr));
```

**Rationale:** Capture full context for next segfault to understand exactly what collides.

---

### 🔍 **INVESTIGATION 2: Reproduce with Stress Test**

Create `tests/test_parallel_stress.rs`:

```rust
#[test]
fn test_parallel_fastga_stress() {
    use std::thread;

    let handles: Vec<_> = (0..10).map(|i| {
        thread::spawn(move || {
            let temp_dir = tempfile::tempdir().unwrap();
            let input = temp_dir.path().join("test.fa.gz");
            fs::copy("data/scerevisiae8.fa.gz", &input).unwrap();

            let result = Command::new("cargo")
                .args(&["run", "--release", "--", input.to_str().unwrap(), "--paf"])
                .output()
                .unwrap();

            assert!(result.status.success(), "Thread {} failed", i);
        })
    }).collect();

    for h in handles {
        h.join().unwrap();
    }
}
```

**Run 100 times:**
```bash
for i in {1..100}; do
    echo "Iteration $i"
    cargo test test_parallel_fastga_stress || break
done
```

**Expected:** Should fail within first 10-20 iterations if race condition is real.

---

### 🔬 **INVESTIGATION 3: Inspect FastGA GDB Creation**

**File:** `deps/fastga/FAtoGDB.c`

**Questions:**
1. Does `FAtoGDB` check for existing GDB files before creating?
2. Is there file locking?
3. What happens if two processes try to create the same GDB?

**Code Review Needed:**
```c
// Look for patterns like:
if (access(gdb_path, F_OK) == 0) {
    // File exists, reuse
} else {
    // Create new - POTENTIAL RACE WINDOW
}
```

**Action:** Manual code review of GDB creation logic.

---

### 🛠️ **FIX 3: Patch FastGA Source (If Needed)**

If GDB creation has race condition:

**Option A: Add File Locking**
```c
int fd = open(gdb_path, O_CREAT | O_EXCL | O_WRONLY, 0644);
if (fd < 0 && errno == EEXIST) {
    // Another process is creating it, wait or error
}
```

**Option B: Use Atomic Operations**
```c
// Use mkstemp for temp file, then rename atomically
```

**Risk:** High (modifying upstream C code)
**Effort:** 2-4 hours (requires C debugging)
**Expected Impact:** 90% chance of fixing if GDB is the issue

---

## Immediate Actions

### ✅ **DONE: Workaround in Place**
- Ubuntu CI runs tests sequentially (`--test-threads=1`)
- Acceptable for now but doubles CI time

### 🔥 **NEXT: Implement FIX 1 + FIX 2**
1. Add unique temp directory via `-P` flag (30 min)
2. Use `tempfile` for .1aln output (15 min)
3. Test on CI with parallel execution restored
4. If passes: Deploy, document, close issue
5. If fails: Proceed to Investigation phase

### 📊 **THEN: Stress Testing**
- Create parallel stress test
- Run locally 100+ times
- Capture any failures with verbose logging

### 🔬 **FINALLY: Upstream Fix (If Needed)**
- Review FastGA C source
- Identify exact race condition
- Patch and submit upstream PR
- Document in UPDATING_FASTGA.md

---

## Testing Checklist

When implementing fixes:

- [ ] Local test with `--test-threads=4` (10 runs)
- [ ] Local test with `--test-threads=8` (10 runs)
- [ ] CI test with parallel execution restored
- [ ] Stress test (100 iterations)
- [ ] Verify no performance regression
- [ ] Verify cleanup works (no leaked temp dirs)
- [ ] Test on macOS (should remain working)
- [ ] Test on Ubuntu (should now work in parallel)

---

## Related Files

**fastga-rs:**
- `src/runner.rs:56` - Temp file creation
- `src/runner.rs:117` - Temp directory flag
- `src/config.rs` - `temp_dir` field
- `.github/workflows/ci.yml` - Sequential execution workaround

**FastGA C Source:**
- `deps/fastga/FastGA.c` - Main entry point
- `deps/fastga/FAtoGDB.c` - FASTA → GDB conversion
- `deps/fastga/align.c` - Multi-threaded alignment

**Sweepga Test Suite:**
- `tests/test_end_to_end.rs` - End-to-end pipeline tests
- `tests/test_large_scale_equivalence.rs` - Large-scale tests
- `tests/test_filtering_modes.rs` - Filtering mode tests

---

## References

- GitHub Actions CI Run (Failed): https://github.com/pangenome/sweepga/actions/runs/18506680661
- GitHub Actions CI Run (Success, Sequential): https://github.com/pangenome/sweepga/actions/runs/18508805953
- FastGA `-P` flag documentation: `FastGA -h` output

---

## Contact

For questions or to report findings:
- Open issue in `fastga-rs` repository
- Tag: `race-condition`, `segfault`, `parallel-testing`

**Last Updated:** October 15, 2025
