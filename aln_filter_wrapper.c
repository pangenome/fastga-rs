#include "aln_filter_wrapper.h"
#include "alncode.h"
#include "GDB.h"
#include <stdlib.h>
#include <string.h>

// Internal handle structure
typedef struct {
    OneFile *of;
    GDB *gdb1;
    GDB *gdb2;
    char *db1_name;
    char *db2_name;
    int64 num_alignments;
    int64 current_idx;
    int tspace;
} AlnHandle;

void* aln_open(const char *path, int64 *num_alignments) {
    AlnHandle *handle = (AlnHandle*)malloc(sizeof(AlnHandle));
    if (!handle) return NULL;

    memset(handle, 0, sizeof(AlnHandle));

    // Open the .1aln file
    char *cpath = NULL;
    handle->of = open_Aln_Read((char*)path, 1,
                               &handle->num_alignments,
                               &handle->tspace,
                               &handle->db1_name,
                               &handle->db2_name,
                               &cpath);

    if (!handle->of) {
        free(handle);
        return NULL;
    }

    // Load the genome databases to get sequence names and lengths
    handle->gdb1 = (GDB*)malloc(sizeof(GDB));
    handle->gdb2 = (GDB*)malloc(sizeof(GDB));

    if (Read_GDB(handle->gdb1, handle->db1_name) != 0) {
        oneFileClose(handle->of);
        free(handle->gdb1);
        free(handle->gdb2);
        free(handle);
        return NULL;
    }

    if (strcmp(handle->db1_name, handle->db2_name) == 0) {
        // Same database for both - just point to the same structure
        handle->gdb2 = handle->gdb1;
    } else {
        if (Read_GDB(handle->gdb2, handle->db2_name) != 0) {
            Close_GDB(handle->gdb1);
            oneFileClose(handle->of);
            free(handle->gdb1);
            free(handle->gdb2);
            free(handle);
            return NULL;
        }
    }

    handle->current_idx = 0;

    if (num_alignments) {
        *num_alignments = handle->num_alignments;
    }

    return handle;
}

int aln_read_record(void *handle_ptr, AlnRecord *record) {
    AlnHandle *handle = (AlnHandle*)handle_ptr;
    if (!handle || !record) return -1;

    if (handle->current_idx >= handle->num_alignments) {
        return -1; // EOF
    }

    // Read the alignment overlap record
    Overlap ovl;
    memset(&ovl, 0, sizeof(Overlap));

    Read_Aln_Overlap(handle->of, &ovl);

    // Skip the trace data - we don't need it for filtering
    Skip_Aln_Trace(handle->of);

    // NOTE: aread/bread in .1aln are CONTIG IDs, not scaffold IDs!
    // We need to map contig -> scaffold to get sequence names and lengths
    // AND adjust coordinates from contig-relative to scaffold-relative

    record->reverse = (ovl.flags & COMP_FLAG) ? 1 : 0;
    record->diffs = ovl.path.diffs;

    // Map query contig ID to scaffold ID and adjust coordinates
    if (ovl.aread >= 0 && ovl.aread < handle->gdb1->ncontig) {
        int ascaff = handle->gdb1->contigs[ovl.aread].scaf;
        int64 aoff = handle->gdb1->contigs[ovl.aread].sbeg;  // Contig offset within scaffold

        record->query_id = ascaff;  // Store scaffold ID
        record->query_len = handle->gdb1->scaffolds[ascaff].slen;
        // Query coordinates are always forward strand
        record->query_start = ovl.path.abpos + aoff;
        record->query_end = ovl.path.aepos + aoff;
    } else {
        record->query_id = -1;
        record->query_len = 0;
        record->query_start = 0;
        record->query_end = 0;
    }

    // Map target contig ID to scaffold ID and adjust coordinates
    if (ovl.bread >= 0 && ovl.bread < handle->gdb2->ncontig) {
        int bscaff = handle->gdb2->contigs[ovl.bread].scaf;
        int64 sbeg = handle->gdb2->contigs[ovl.bread].sbeg;
        int64 clen = handle->gdb2->contigs[ovl.bread].clen;

        record->target_id = bscaff;  // Store scaffold ID
        record->target_len = handle->gdb2->scaffolds[bscaff].slen;

        // Target coordinates depend on strand
        if (record->reverse) {
            // Reverse strand: coordinates relative to reverse complement
            int64 boff = sbeg + clen;
            record->target_start = boff - ovl.path.bepos;
            record->target_end = boff - ovl.path.bbpos;
        } else {
            // Forward strand: add contig offset
            record->target_start = ovl.path.bbpos + sbeg;
            record->target_end = ovl.path.bepos + sbeg;
        }
    } else {
        record->target_id = -1;
        record->target_len = 0;
        record->target_start = 0;
        record->target_end = 0;
    }

    handle->current_idx++;

    return 0;
}

const char* aln_get_seq_name(void *handle_ptr, int64 seq_id, int which_db) {
    AlnHandle *handle = (AlnHandle*)handle_ptr;
    if (!handle) return NULL;

    GDB *gdb = (which_db == 0) ? handle->gdb1 : handle->gdb2;

    if (seq_id < 0 || seq_id >= gdb->nscaff) {
        return NULL;
    }

    // Get the scaffold header from the headers block
    int64 hoff = gdb->scaffolds[seq_id].hoff;
    if (hoff < 0 || hoff >= gdb->hdrtot) {
        return NULL;
    }

    return gdb->headers + hoff;
}

void aln_close(void *handle_ptr) {
    AlnHandle *handle = (AlnHandle*)handle_ptr;
    if (!handle) return;

    if (handle->of) {
        oneFileClose(handle->of);
    }

    if (handle->gdb1) {
        Close_GDB(handle->gdb1);
        free(handle->gdb1);
    }

    if (handle->gdb2 && handle->gdb2 != handle->gdb1) {
        Close_GDB(handle->gdb2);
        free(handle->gdb2);
    }

    free(handle);
}
