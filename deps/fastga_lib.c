/* FastGA Library Interface for Streaming Alignments
 *
 * This file provides a library interface to FastGA that allows intercepting
 * alignments as they are generated, enabling streaming processing without
 * writing intermediate files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "GDB.h"
#include "align.h"
#include "alncode.h"
#include "ONElib.h"

// Callback function type for streaming alignments
typedef int (*alignment_callback)(
    void* user_data,
    const char* query_name,
    int query_len,
    int query_start,
    int query_end,
    const char* target_name,
    int target_len,
    int target_start,
    int target_end,
    int strand,  // 0 for forward, 1 for reverse
    const char* cigar,
    int matches,
    int mismatches,
    int gaps
);

// Context for streaming alignments
typedef struct {
    alignment_callback callback;
    void* user_data;
    int alignment_count;
    int filtered_count;
} StreamContext;

// Global streaming context (will be thread-local in production)
static StreamContext* g_stream_context = NULL;

// Initialize streaming context
StreamContext* fastga_stream_init(alignment_callback callback, void* user_data) {
    StreamContext* ctx = (StreamContext*)malloc(sizeof(StreamContext));
    if (ctx) {
        ctx->callback = callback;
        ctx->user_data = user_data;
        ctx->alignment_count = 0;
        ctx->filtered_count = 0;
    }
    return ctx;
}

// Free streaming context
void fastga_stream_free(StreamContext* ctx) {
    if (ctx) {
        free(ctx);
    }
}

// Convert alignment to CIGAR with extended operators
char* alignment_to_extended_cigar(Alignment* aln) {
    static char cigar_buffer[65536];  // Static buffer for CIGAR string
    cigar_buffer[0] = '\0';

    Path* path = aln->path;
    if (!path || !path->trace) {
        return cigar_buffer;
    }

    // Simplified CIGAR generation - in production this would use the full
    // trace point decoding logic from ALNtoPAF.c
    int trace_len = path->tlen;
    int32_t* trace = (int32_t*)path->trace;

    char* ptr = cigar_buffer;
    int remaining = sizeof(cigar_buffer) - 1;

    // Generate extended CIGAR based on trace points
    int matches = 0, mismatches = 0, insertions = 0, deletions = 0;

    for (int i = 0; i < trace_len && remaining > 20; i++) {
        int32_t val = trace[i];
        if (val < 0) {
            // Deletion
            deletions++;
            if (matches > 0) {
                int written = snprintf(ptr, remaining, "%d=", matches);
                ptr += written;
                remaining -= written;
                matches = 0;
            }
            if (mismatches > 0) {
                int written = snprintf(ptr, remaining, "%dX", mismatches);
                ptr += written;
                remaining -= written;
                mismatches = 0;
            }
        } else {
            // Match or mismatch - simplified for now
            matches++;
            if (deletions > 0) {
                int written = snprintf(ptr, remaining, "%dD", deletions);
                ptr += written;
                remaining -= written;
                deletions = 0;
            }
            if (insertions > 0) {
                int written = snprintf(ptr, remaining, "%dI", insertions);
                ptr += written;
                remaining -= written;
                insertions = 0;
            }
        }
    }

    // Flush remaining operations
    if (matches > 0) {
        snprintf(ptr, remaining, "%d=", matches);
    } else if (mismatches > 0) {
        snprintf(ptr, remaining, "%dX", mismatches);
    } else if (deletions > 0) {
        snprintf(ptr, remaining, "%dD", deletions);
    } else if (insertions > 0) {
        snprintf(ptr, remaining, "%dI", insertions);
    }

    return cigar_buffer;
}

// Hook function to intercept alignments
int fastga_process_alignment(GDB* gdb1, GDB* gdb2, Overlap* ovl, Alignment* aln) {
    if (!g_stream_context || !g_stream_context->callback) {
        return 1;  // No callback, continue normally
    }

    g_stream_context->alignment_count++;

    // Extract alignment information
    char* query_name = gdb1->scafolds[ovl->aread].name;
    char* target_name = gdb2->scafolds[ovl->bread].name;

    int query_len = gdb1->scafolds[ovl->aread].clen;
    int target_len = gdb2->scafolds[ovl->bread].clen;

    Path* path = &ovl->path;
    int query_start = path->abpos;
    int query_end = path->aepos;
    int target_start = path->bbpos;
    int target_end = path->bepos;

    int strand = (ovl->flags & 0x1) ? 1 : 0;  // COMP_FLAG

    // Generate extended CIGAR
    char* cigar = alignment_to_extended_cigar(aln);

    // Calculate statistics
    int matches = 0, mismatches = 0, gaps = 0;
    // These would be calculated from the actual alignment
    matches = (query_end - query_start) * 0.9;  // Simplified
    mismatches = (query_end - query_start) * 0.1;

    // Call user callback
    int result = g_stream_context->callback(
        g_stream_context->user_data,
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
        gaps
    );

    if (result == 0) {
        g_stream_context->filtered_count++;
        return 0;  // Skip this alignment
    }

    return 1;  // Continue processing
}

// Main entry point for streaming alignment
int fastga_align_streaming(
    const char* genome1_path,
    const char* genome2_path,
    alignment_callback callback,
    void* user_data,
    int num_threads,
    int min_length,
    double min_identity
) {
    // Set up streaming context
    StreamContext* ctx = fastga_stream_init(callback, user_data);
    if (!ctx) {
        return -1;
    }

    g_stream_context = ctx;

    // Here we would call the actual FastGA alignment functions
    // For now, this is a placeholder that demonstrates the interface

    // In production, this would:
    // 1. Load the GDB files
    // 2. Build or load indices
    // 3. Run the alignment with our callback hook
    // 4. Clean up resources

    int result = 0;

    // Report statistics
    fprintf(stderr, "Processed %d alignments, kept %d\n",
            ctx->alignment_count,
            ctx->alignment_count - ctx->filtered_count);

    // Clean up
    g_stream_context = NULL;
    fastga_stream_free(ctx);

    return result;
}