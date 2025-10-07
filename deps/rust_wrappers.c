// Wrapper functions for FastGA utilities that can be called from Rust
// These are simple stubs that just return -1 to indicate they should use the main functions

#include <stdio.h>

// Wrapper for FAtoGDB - not implemented, use fatogdb_main instead
int fatogdb_wrapper(const char *input_path) {
    fprintf(stderr, "[Wrapper] FAtoGDB wrapper not implemented, use fatogdb_main\n");
    return -1;
}

// Wrapper for GIXmake
// For now, we'll have to use the full program since GIXmake is complex
int gixmake_wrapper(const char *gdb_path, int threads, const char *temp_dir, int freq) {
    fprintf(stderr, "[Wrapper] GIXmake wrapper called for %s\n", gdb_path);
    // This would require extracting the core GIXmake logic
    // For now, return error to indicate we need to use the binary
    return -1;
}