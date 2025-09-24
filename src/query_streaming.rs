/// Efficient query streaming using FastGA's natural processing order
///
/// Key insight: FastGA already processes queries sequentially!
/// It completes all alignments for query[0] before starting query[1]

use std::path::Path;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead};
use anyhow::Result;
use crate::{Alignment, Config};

pub struct StreamingQueryAligner {
    target_gdb: String,
    config: Config,
}

impl StreamingQueryAligner {
    pub fn new(target_db: impl AsRef<Path>, config: Config) -> Result<Self> {
        let target_gdb = crate::ffi::prepare_database(target_db.as_ref())?;
        Ok(Self { target_gdb, config })
    }

    /// Stream alignments with query boundaries guaranteed by FastGA's processing order
    pub fn align_queries_streaming<F>(
        &self,
        queries_fasta: impl AsRef<Path>,
        mut on_query_complete: F,
    ) -> Result<()>
    where
        F: FnMut(String, Vec<Alignment>) -> bool,  // (query_name, alignments) -> continue?
    {
        // Convert queries to GDB once
        let query_gdb = crate::ffi::prepare_database(queries_fasta.as_ref())?;

        // Run FastGA with streaming output
        let mut child = Command::new("FastGA")
            .arg("-T").arg(self.config.num_threads.to_string())
            .arg("-pafx")  // PAF with extended CIGAR
            .arg(&query_gdb)
            .arg(&self.target_gdb)
            .stdout(Stdio::piped())
            .spawn()?;

        let stdout = child.stdout.take().expect("Failed to capture stdout");
        let reader = BufReader::new(stdout);

        let mut current_query: Option<String> = None;
        let mut current_alignments: Vec<Alignment> = Vec::new();

        // Parse PAF output line by line
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Parse PAF line to get query name
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                continue;  // Invalid PAF line
            }

            let query_name = fields[0].to_string();

            // Check if we've moved to a new query
            if let Some(ref prev_query) = current_query {
                if query_name != *prev_query {
                    // Previous query is complete!
                    if !on_query_complete(prev_query.clone(), current_alignments.clone()) {
                        // User wants to stop
                        child.kill()?;
                        return Ok(());
                    }
                    current_alignments.clear();
                }
            }

            current_query = Some(query_name);

            // Parse alignment and add to current batch
            let alignment = parse_paf_line(&line)?;
            current_alignments.push(alignment);
        }

        // Handle last query
        if let Some(query) = current_query {
            on_query_complete(query, current_alignments);
        }

        child.wait()?;
        Ok(())
    }

    /// Even simpler: Just track when each query is done
    pub fn process_with_boundaries<F, G>(
        &self,
        queries_fasta: impl AsRef<Path>,
        mut on_alignment: F,
        mut on_query_done: G,
    ) -> Result<()>
    where
        F: FnMut(Alignment),
        G: FnMut(String, usize),  // (query_name, alignment_count)
    {
        let query_gdb = crate::ffi::prepare_database(queries_fasta.as_ref())?;

        let output = Command::new("FastGA")
            .arg("-T").arg(self.config.num_threads.to_string())
            .arg("-pafx")
            .arg(&query_gdb)
            .arg(&self.target_gdb)
            .output()?;

        let paf = String::from_utf8(output.stdout)?;

        let mut current_query: Option<String> = None;
        let mut count = 0;

        for line in paf.lines() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let alignment = parse_paf_line(line)?;

            // Detect query boundary
            if let Some(ref prev) = current_query {
                if alignment.query_name != *prev {
                    on_query_done(prev.clone(), count);
                    count = 0;
                }
            }

            current_query = Some(alignment.query_name.clone());
            count += 1;
            on_alignment(alignment);
        }

        // Last query
        if let Some(query) = current_query {
            on_query_done(query, count);
        }

        Ok(())
    }
}

fn parse_paf_line(line: &str) -> Result<Alignment> {
    let fields: Vec<&str> = line.split('\t').collect();

    // PAF format:
    // 0: query_name
    // 1: query_len
    // 2: query_start
    // 3: query_end
    // 4: strand
    // 5: target_name
    // 6: target_len
    // 7: target_start
    // 8: target_end
    // 9: num_matches
    // 10: block_len
    // 11: mapq
    // 12+: tags (including cg:Z:cigar)

    let mut alignment = Alignment {
        query_name: fields[0].to_string(),
        query_len: fields[1].parse()?,
        query_start: fields[2].parse()?,
        query_end: fields[3].parse()?,
        strand: fields[4].chars().next().unwrap_or('+'),
        target_name: fields[5].to_string(),
        target_len: fields[6].parse()?,
        target_start: fields[7].parse()?,
        target_end: fields[8].parse()?,
        residue_matches: fields[9].parse()?,
        block_length: fields[10].parse()?,
        mapping_quality: fields[11].parse()?,
        cigar: String::new(),
        ..Default::default()
    };

    // Parse tags for CIGAR
    for tag in &fields[12..] {
        if tag.starts_with("cg:Z:") {
            alignment.cigar = tag[5..].to_string();
            break;
        }
    }

    Ok(alignment)
}

/// Example usage:
/// ```
/// let aligner = StreamingQueryAligner::new("target.fa", Config::default())?;
///
/// // Process queries with guaranteed completeness
/// aligner.align_queries_streaming("queries.fa", |query_name, alignments| {
///     println!("Query {} complete with {} alignments", query_name, alignments.len());
///     // Process this query's results
///     true  // Continue to next query
/// })?;
/// ```