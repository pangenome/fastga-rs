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
                               NULL,  // gdb1 - will be loaded later
                               NULL,  // gdb2 - will be loaded later
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

// === Writer implementation ===

typedef struct {
    OneFile *of;
    GDB *gdb1;
    GDB *gdb2;
} AlnWriter;

void* aln_create(const char *path, const char *gdb1_path, const char *gdb2_path) {
    AlnWriter *writer = (AlnWriter*)malloc(sizeof(AlnWriter));
    if (!writer) return NULL;

    writer->gdb1 = NULL;
    writer->gdb2 = NULL;
    writer->of = NULL;

    // Load GDB files for sequence metadata
    writer->gdb1 = (GDB*)malloc(sizeof(GDB));
    if (!writer->gdb1) goto error;

    if (Read_GDB(writer->gdb1, (char*)gdb1_path) < 0) {
        fprintf(stderr, "Failed to read GDB: %s\n", gdb1_path);
        goto error;
    }

    // Check if same file for both databases
    if (strcmp(gdb1_path, gdb2_path) == 0) {
        writer->gdb2 = writer->gdb1;
    } else {
        writer->gdb2 = (GDB*)malloc(sizeof(GDB));
        if (!writer->gdb2) goto error;

        if (Read_GDB(writer->gdb2, (char*)gdb2_path) < 0) {
            fprintf(stderr, "Failed to read GDB: %s\n", gdb2_path);
            goto error;
        }
    }

    // Create ONEcode schema and file
    OneSchema *schema = make_Aln_Schema();
    if (!schema) {
        fprintf(stderr, "Failed to create .1aln schema\n");
        goto error;
    }

    writer->of = oneFileOpenWriteNew((char*)path, schema, "aln", 1, 1);
    if (!writer->of) {
        fprintf(stderr, "Failed to create .1aln file: %s\n", path);
        oneSchemaDestroy(schema);
        goto error;
    }

    // Write provenance
    oneAddProvenance(writer->of, "sweepga", "0.1.0", "sweepga filter");

    // Write references to GDB files
    oneAddReference(writer->of, (char*)gdb1_path, 1);
    oneAddReference(writer->of, (char*)gdb2_path, 2);

    // Write trace point spacing (required by schema)
    // Use 100 as default (same as FastGA)
    oneInt(writer->of, 0) = 100;
    oneWriteLine(writer->of, 't', 0, 0);

    return writer;

error:
    if (writer) {
        if (writer->gdb1) {
            Close_GDB(writer->gdb1);
            free(writer->gdb1);
        }
        if (writer->gdb2 && writer->gdb2 != writer->gdb1) {
            Close_GDB(writer->gdb2);
            free(writer->gdb2);
        }
        if (writer->of) {
            oneFileClose(writer->of);
        }
        free(writer);
    }
    return NULL;
}

int aln_write_record(void *handle_ptr, const AlnRecord *rec) {
    AlnWriter *writer = (AlnWriter*)handle_ptr;
    if (!writer || !writer->of || !rec) return -1;

    // Map scaffold IDs back to contig IDs
    // For now, assume scaffold ID == contig ID (simplified)
    // A full implementation would need reverse lookup

    int64 aread = rec->query_id;
    int64 bread = rec->target_id;

    // Write alignment record: O A 6 3 INT (aread, abpos, aepos, bread, bbpos, bepos)
    oneInt(writer->of, 0) = aread;
    oneInt(writer->of, 1) = rec->query_start;
    oneInt(writer->of, 2) = rec->query_end;
    oneInt(writer->of, 3) = bread;
    oneInt(writer->of, 4) = rec->target_start;
    oneInt(writer->of, 5) = rec->target_end;
    oneWriteLine(writer->of, 'A', 0, 0);

    // Write reverse flag if set: D R 0
    if (rec->reverse) {
        oneWriteLine(writer->of, 'R', 0, 0);
    }

    // Write differences: D D 1 3 INT
    oneInt(writer->of, 0) = rec->diffs;
    oneWriteLine(writer->of, 'D', 0, 0);

    // Write sequence lengths: D L 2 3 INT 3 INT
    oneInt(writer->of, 0) = rec->query_len;
    oneInt(writer->of, 1) = rec->target_len;
    oneWriteLine(writer->of, 'L', 0, 0);

    // Write trace points for ALNtoPAF compatibility
    // T line: trace points in target (b) sequence
    // We write a single trace point at the end of the alignment
    int64 tlen = rec->target_end - rec->target_start;
    int64 trace_points[1] = { tlen };
    oneWriteLine(writer->of, 'T', 1, trace_points);

    // X line: number of diffs per trace interval
    // All diffs are in the single interval
    int64 diffs_list[1] = { rec->diffs };
    oneWriteLine(writer->of, 'X', 1, diffs_list);

    return 0;
}

void aln_close_writer(void *handle_ptr) {
    AlnWriter *writer = (AlnWriter*)handle_ptr;
    if (!writer) return;

    if (writer->of) {
        oneFileClose(writer->of);
    }

    if (writer->gdb1) {
        Close_GDB(writer->gdb1);
        free(writer->gdb1);
    }

    if (writer->gdb2 && writer->gdb2 != writer->gdb1) {
        Close_GDB(writer->gdb2);
        free(writer->gdb2);
    }

    free(writer);
}
