// FastGA Query-Specific API Extension
// Allows aligning one query at a time against a database

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "GDB.h"
#include "align.h"
#include "gene_core.h"

// External functions from FastGA.c we need
extern void align_contigs(uint8_t *beg, uint8_t *end, int swide, int ctg1, int ctg2,
                         void *pair);  // Contig_Bundle

// Structure to hold alignment results for a single query
typedef struct {
    int query_id;
    int target_id;
    int query_start;
    int query_end;
    int target_start;
    int target_end;
    double identity;
    char *cigar;
    int strand;  // 0 = forward, 1 = reverse
} QueryAlignment;

typedef struct {
    QueryAlignment *alignments;
    int count;
    int capacity;
} QueryAlignmentSet;

// Callback function type for streaming alignments
typedef int (*AlignmentCallback)(const QueryAlignment *aln, void *user_data);

// Main API function: Align a single query sequence against a target database
// Returns all alignments for this query
QueryAlignmentSet* align_single_query(
    GDB *query_gdb,       // GDB containing single query sequence
    int query_idx,        // Index of query in the GDB (usually 0)
    GDB *target_gdb,      // Target database
    Work_Data *work,      // Alignment work data
    Align_Spec *spec,     // Alignment specification
    AlignmentCallback callback,  // Optional callback for streaming
    void *user_data       // User data for callback
) {
    QueryAlignmentSet *result = (QueryAlignmentSet*)calloc(1, sizeof(QueryAlignmentSet));
    if (!result) return NULL;

    result->capacity = 100;
    result->alignments = (QueryAlignment*)malloc(result->capacity * sizeof(QueryAlignment));
    if (!result->alignments) {
        free(result);
        return NULL;
    }

    // Get query contig
    if (query_idx >= query_gdb->ncontig) {
        free(result->alignments);
        free(result);
        return NULL;
    }

    Contig *query_contig = &query_gdb->contigs[query_idx];

    // Align query against all targets
    for (int target_idx = 0; target_idx < target_gdb->ncontig; target_idx++) {
        Contig *target_contig = &target_gdb->contigs[target_idx];

        // Skip if either is masked/invalid
        if (query_contig->boff < 0 || target_contig->boff < 0)
            continue;

        // TODO: Call FastGA's alignment machinery here
        // This would involve setting up the kmer matching and calling align_contigs
        // For now, this is a placeholder showing the API structure

        // If alignment found, add to results
        QueryAlignment aln = {0};
        aln.query_id = query_idx;
        aln.target_id = target_idx;
        // ... fill in alignment details ...

        // Call callback if provided
        if (callback) {
            if (!callback(&aln, user_data)) {
                // Callback returned 0, stop processing
                break;
            }
        }

        // Add to result set
        if (result->count >= result->capacity) {
            result->capacity *= 2;
            QueryAlignment *new_alns = (QueryAlignment*)realloc(
                result->alignments,
                result->capacity * sizeof(QueryAlignment)
            );
            if (!new_alns) {
                free(result->alignments);
                free(result);
                return NULL;
            }
            result->alignments = new_alns;
        }
        result->alignments[result->count++] = aln;
    }

    return result;
}

// Free alignment results
void free_query_alignments(QueryAlignmentSet *set) {
    if (set) {
        if (set->alignments) {
            for (int i = 0; i < set->count; i++) {
                if (set->alignments[i].cigar)
                    free(set->alignments[i].cigar);
            }
            free(set->alignments);
        }
        free(set);
    }
}

// Batch API: Process multiple queries with guaranteed completeness
// This ensures each query is fully processed before moving to the next
int align_queries_streaming(
    GDB *query_gdb,       // GDB with multiple query sequences
    GDB *target_gdb,      // Target database
    Work_Data *work,      // Alignment work data
    Align_Spec *spec,     // Alignment specification
    AlignmentCallback callback,  // Callback for each alignment
    void *user_data       // User data for callback
) {
    int total_alignments = 0;

    // Process each query in order
    for (int q = 0; q < query_gdb->ncontig; q++) {
        QueryAlignmentSet *results = align_single_query(
            query_gdb, q, target_gdb,
            work, spec, callback, user_data
        );

        if (results) {
            total_alignments += results->count;
            free_query_alignments(results);
        }
    }

    return total_alignments;
}