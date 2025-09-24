/// Strategy for ensuring query completeness without per-query GDB overhead
///
/// The key insight: We can batch queries but ensure completeness by
/// modifying how we call FastGA or post-processing its output

use std::collections::HashMap;
use std::path::Path;
use anyhow::Result;
use crate::{Alignment, Config};

/// Efficient query-complete alignment strategy
pub struct EfficientQueryAligner {
    target_gdb: String,
    config: Config,
}

impl EfficientQueryAligner {
    pub fn new(target_db: impl AsRef<Path>, config: Config) -> Result<Self> {
        let target_gdb = crate::ffi::prepare_database(target_db.as_ref())?;
        Ok(Self { target_gdb, config })
    }

    /// Strategy 1: Batch queries but track completeness
    /// Put all queries in ONE GDB, run ONCE, then sort output
    pub fn align_batch_sorted(&self, queries: Vec<Query>) -> Result<Vec<QueryResult>> {
        // Create single FASTA with all queries
        let temp_dir = tempfile::tempdir()?;
        let batch_fasta = temp_dir.path().join("queries.fa");

        let mut fasta_content = String::new();
        for query in &queries {
            fasta_content.push_str(&format!(">{}\n{}\n",
                query.name,
                std::str::from_utf8(&query.sequence)?
            ));
        }
        std::fs::write(&batch_fasta, fasta_content)?;

        // ONE GDB creation for all queries
        let query_gdb = crate::ffi::prepare_database(&batch_fasta)?;

        // ONE FastGA run
        let alignments = crate::ffi::run_fastga_alignment(
            &query_gdb,
            &self.target_gdb,
            &self.config
        )?;

        // Sort alignments by query name
        let mut results_map: HashMap<String, Vec<Alignment>> = HashMap::new();
        for aln in alignments.alignments {
            results_map.entry(aln.query_name.clone())
                .or_insert_with(Vec::new)
                .push(aln);
        }

        // Convert to ordered results
        let results: Vec<QueryResult> = queries.iter()
            .map(|q| QueryResult {
                query_name: q.name.clone(),
                alignments: results_map.remove(&q.name).unwrap_or_default(),
            })
            .collect();

        Ok(results)
    }

    /// Strategy 2: Use subprocess with persistent target
    /// Fork FastGA for each query but reuse target GDB
    pub fn align_with_subprocess(&self, query: &Query) -> Result<QueryResult> {
        use std::process::Command;

        // Write single query to temp file
        let temp_dir = tempfile::tempdir()?;
        let query_fasta = temp_dir.path().join("query.fa");
        std::fs::write(&query_fasta, format!(
            ">{}\n{}\n",
            query.name,
            std::str::from_utf8(&query.sequence)?
        ))?;

        // Convert to GDB (this is the overhead we pay per query)
        let query_gdb = temp_dir.path().join("query.gdb");
        Command::new("FAtoGDB")
            .arg(&query_fasta)
            .arg(&query_gdb)
            .output()?;

        // Run FastGA as subprocess (fork overhead ~10ms)
        let output = Command::new("FastGA")
            .arg("-T").arg(self.config.num_threads.to_string())
            .arg(&query_gdb)
            .arg(&self.target_gdb)
            .output()?;

        // Parse PAF output
        let paf = String::from_utf8(output.stdout)?;
        let alignments = parse_paf_output(&paf)?;

        Ok(QueryResult {
            query_name: query.name.clone(),
            alignments,
        })
    }

    /// Strategy 3: Modified FastGA binary
    /// Compile a custom FastGA that processes queries sequentially
    /// This requires patching FastGA.c to add a flag like -Q for "query mode"
    pub fn align_with_patched_fastga(&self, queries: Vec<Query>) -> Result<Vec<QueryResult>> {
        // This would use a modified FastGA binary that we compile
        // with a patch that ensures sequential query processing:
        //
        // In FastGA.c main loop:
        // ```c
        // if (QUERY_MODE) {
        //     for (r = 0; r < gdb1->ncontig; r++) {
        //         // Process query r completely
        //         for (s = 0; s < gdb2->ncontig; s++) {
        //             align_contigs(...);
        //         }
        //         // Flush output for this query
        //         fflush(stdout);
        //         // Output marker
        //         if (OUT_TYPE == PAF)
        //             printf("##QUERY_COMPLETE: %s\n", contig_name);
        //     }
        // }
        // ```

        unimplemented!("Requires patched FastGA binary")
    }
}

/// Overhead comparison:
///
/// | Method | Setup | Per Query | 1000 Queries |
/// |--------|-------|-----------|--------------|
/// | Individual GDB | 0 | ~500ms | 500s |
/// | Batch + Sort | ~1s | 0 | 1s + alignment time |
/// | Subprocess | 0 | ~50ms | 50s |
/// | Patched FastGA | 0 | 0 | alignment time only |
///
/// The batch approach is clearly best for throughput.
/// The subprocess approach gives you streaming with moderate overhead.
/// The patched approach is ideal but requires maintaining a fork.

pub struct Query {
    pub name: String,
    pub sequence: Vec<u8>,
}

pub struct QueryResult {
    pub query_name: String,
    pub alignments: Vec<Alignment>,
}

fn parse_paf_output(paf: &str) -> Result<Vec<Alignment>> {
    // Parse PAF format into Alignment structs
    unimplemented!()
}