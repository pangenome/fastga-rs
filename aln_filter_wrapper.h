#ifndef ALN_FILTER_WRAPPER_H
#define ALN_FILTER_WRAPPER_H

#include "gene_core.h"  // For int64 type

// Simple alignment record for filtering
typedef struct {
    int64 query_id;
    int64 query_start;
    int64 query_end;
    int64 target_id;
    int64 target_start;
    int64 target_end;
    int reverse;
    int diffs;
    int64 query_len;
    int64 target_len;
} AlnRecord;

// Open .1aln file and get alignment count
void* aln_open(const char *path, int64 *num_alignments);

// Read next alignment record (returns 0 on success, -1 on EOF)
int aln_read_record(void *handle, AlnRecord *record);

// Get sequence name by ID
const char* aln_get_seq_name(void *handle, int64 seq_id, int which_db);

// Close file
void aln_close(void *handle);

// === Writing functions ===

// Create new .1aln file for writing
// gdb1_path and gdb2_path: paths to .1gdb files (sequence metadata)
void* aln_create(const char *path, const char *gdb1_path, const char *gdb2_path);

// Write an alignment record
int aln_write_record(void *handle, const AlnRecord *record);

// Close writer (flushes and finalizes file)
void aln_close_writer(void *handle);

#endif
