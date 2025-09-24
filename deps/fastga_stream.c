/* FastGA Streaming Integration
 *
 * This file modifies FastGA to expose alignment generation as a library
 * with callback support for streaming processing.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GDB.h"
#include "align.h"
#include "alncode.h"

// Global callback for streaming alignments
typedef int (*AlignmentCallback)(
    void* user_data,
    const Overlap* ovl,
    const Alignment* aln,
    const char* query_name,
    const char* target_name
);

static AlignmentCallback g_callback = NULL;
static void* g_user_data = NULL;

// Modified version of FastGA's main alignment loop
// This is extracted from FastGA.c and modified to support callbacks
int fastga_align_with_callback(
    const char* genome1_path,
    const char* genome2_path,
    AlignmentCallback callback,
    void* user_data,
    int num_threads,
    int min_length,
    double min_identity
) {
    // Store callback globally
    g_callback = callback;
    g_user_data = user_data;

    // This would contain the core FastGA alignment logic
    // For now, this is a simplified version

    GDB gdb1, gdb2;

    // Load genome databases
    if (Read_GDB(&gdb1, genome1_path) != 0) {
        fprintf(stderr, "Failed to read genome1\n");
        return -1;
    }

    if (Read_GDB(&gdb2, genome2_path) != 0) {
        fprintf(stderr, "Failed to read genome2\n");
        Close_GDB(&gdb1);
        return -1;
    }

    // TODO: Perform actual alignment using FastGA's algorithm
    // This would involve:
    // 1. Building/loading indices
    // 2. Finding seeds
    // 3. Extending alignments
    // 4. Calling our callback for each alignment

    Close_GDB(&gdb2);
    Close_GDB(&gdb1);

    return 0;
}

// Hook function to intercept alignment output
// This replaces Write_Aln_Overlap in the streaming version
int Stream_Aln_Overlap(OneFile* of, Overlap* ovl, Alignment* aln, GDB* gdb1, GDB* gdb2) {
    if (g_callback && ovl && aln) {
        // Get sequence names
        const char* query_name = gdb1->scaffolds[ovl->aread].name;
        const char* target_name = gdb2->scaffolds[ovl->bread].name;

        // Call user callback
        int result = g_callback(g_user_data, ovl, aln, query_name, target_name);

        // If callback returns 0, skip this alignment
        if (result == 0) {
            return 0;
        }
    }

    // If no callback or callback returned 1, write normally
    if (of) {
        Write_Aln_Overlap(of, ovl);
    }

    return 1;
}