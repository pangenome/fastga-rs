# FastGA C Integration Plan

## Current Situation
We have a working Rust wrapper that uses FastGA as a subprocess. This works but doesn't give us the fine-grained control you want for streaming alignments.

## Why True C Integration is Complex

### 1. FastGA's Architecture
- FastGA is written as a standalone program, not a library
- Heavy use of global variables (NTHREADS, ALIGN_MIN, etc.)
- Main function does extensive setup and coordination
- Outputs to files, not to callbacks

### 2. The Alignment Pipeline
```
FastGA.c:main()
  ├── Load/create GDB files
  ├── Build/load GIX indices
  ├── Find seeds (compare_thread)
  ├── Chain seeds
  ├── Generate alignments
  └── Write to .1aln file → Write_Aln_Overlap() at line 4096
```

### 3. Integration Points
The key function where alignments are output is:
```c
// FastGA.c:4096
Write_Aln_Overlap(of, ov);
Write_Aln_Trace(of, src->ptr, tsize, trace64);
```

## Options for True C Integration

### Option 1: Minimal Modification (Recommended)
1. Copy FastGA.c to FastGA_lib.c
2. Replace main() with library entry point
3. Add callback at Write_Aln_Overlap point
4. Compile as static library
5. Link with Rust via FFI

**Pros:**
- Minimal changes to FastGA code
- Preserves all FastGA functionality
- Can intercept alignments as generated

**Cons:**
- Still uses global variables
- Not thread-safe for multiple concurrent alignments

### Option 2: Full Refactor
1. Refactor FastGA to remove globals
2. Create proper library API
3. Make thread-safe
4. Add streaming callbacks throughout

**Pros:**
- Clean library design
- Thread-safe
- Maximum flexibility

**Cons:**
- Massive undertaking
- Risk of introducing bugs
- Diverges from upstream FastGA

### Option 3: Named Pipe/Socket IPC
1. Modify FastGA to write to named pipe
2. Rust reads from pipe in real-time
3. Process alignments as they arrive

**Pros:**
- Minimal FastGA changes
- Process isolation
- Can work across network

**Cons:**
- IPC overhead
- Platform-specific code
- Still not in-process

## Recommended Approach

### Phase 1: Working Prototype (Current)
✅ Use subprocess with PAF parsing (DONE)

### Phase 2: Minimal C Integration
1. Create FastGA_lib.c with:
   ```c
   typedef int (*AlignCallback)(void*, Overlap*, const char*, const char*);

   int fastga_lib_align(
       const char* genome1,
       const char* genome2,
       AlignCallback callback,
       void* user_data,
       FastGAConfig* config
   );
   ```

2. Hook into Write_Aln_Overlap:
   ```c
   if (g_callback) {
       int keep = g_callback(user_data, ov, query_name, target_name);
       if (!keep) continue;  // Skip alignment
   }
   Write_Aln_Overlap(of, ov);
   ```

3. Compile as library and link with Rust

### Phase 3: Production Integration
- Work with FastGA maintainers to add official library support
- Contribute patches upstream
- Maintain fork if necessary

## Implementation Steps

1. **Extract Core Logic** (2-3 days)
   - Copy FastGA.c → FastGA_lib.c
   - Extract main() logic into library function
   - Add callback mechanism

2. **Build System** (1 day)
   - Modify build.rs to compile FastGA_lib.c
   - Link as static library
   - Handle dependencies (pthread, zlib)

3. **FFI Bindings** (1-2 days)
   - Define Overlap struct in Rust
   - Create safe wrappers
   - Handle callback lifetime issues

4. **Testing** (2-3 days)
   - Verify alignment output matches original
   - Test callback filtering
   - Benchmark performance

5. **Documentation** (1 day)
   - Document C modifications
   - API usage examples
   - Performance characteristics

## Complexity Assessment

**Estimated effort:** 1-2 weeks for production-ready integration

**Technical challenges:**
- Global state management
- Memory ownership across FFI boundary
- Error handling in callbacks
- Alignment data structure complexity

**Maintenance burden:**
- Need to track upstream FastGA changes
- Potential merge conflicts
- Testing requirements double

## Recommendation

Given the complexity, I recommend:

1. **Short term:** Keep current subprocess approach for most use cases
2. **Medium term:** Implement minimal C integration for streaming
3. **Long term:** Work with FastGA maintainers on official library API

The subprocess approach with PAF parsing works and is maintainable. True C integration should only be pursued if the performance benefits justify the complexity.