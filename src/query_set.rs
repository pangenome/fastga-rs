use crate::{Alignment, Config};
use anyhow::Result;
use std::io::{BufRead, BufReader};
use std::path::Path;
/// Query-centric alignment API
///
/// Provides guaranteed complete alignment sets for each query sequence,
/// with automatic backpressure to prevent memory overflow.
use std::process::{Command, Stdio};
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread;

/// A complete set of alignments for a single query sequence.
///
/// This represents ALL alignments found for one query against the entire target database.
/// FastGA guarantees that all alignments for a query are reported before moving to the next query.
#[derive(Debug, Clone)]
pub struct QueryAlignmentSet {
    /// Name of the query sequence
    pub query_name: String,

    /// All alignments found for this query
    pub alignments: Vec<Alignment>,
}

impl QueryAlignmentSet {
    /// Get the number of alignments found for this query
    pub fn alignment_count(&self) -> usize {
        self.alignments.len()
    }

    /// Check if this query had any alignments
    pub fn has_alignments(&self) -> bool {
        !self.alignments.is_empty()
    }

    /// Get the best alignment by identity score
    pub fn best_by_identity(&self) -> Option<&Alignment> {
        self.alignments
            .iter()
            .max_by(|a, b| a.identity().partial_cmp(&b.identity()).unwrap())
    }

    /// Filter alignments by minimum identity threshold
    pub fn filter_by_identity(&self, min_identity: f64) -> Vec<&Alignment> {
        self.alignments
            .iter()
            .filter(|a| a.identity() >= min_identity)
            .collect()
    }

    /// Get total query coverage (non-overlapping bases)
    pub fn query_coverage(&self) -> usize {
        // Simple approximation - would need interval merging for exact coverage
        self.alignments
            .iter()
            .map(|a| a.query_end - a.query_start)
            .sum()
    }

    /// Get all unique target sequences hit
    pub fn target_names(&self) -> Vec<String> {
        let mut targets: Vec<_> = self
            .alignments
            .iter()
            .map(|a| a.target_name.clone())
            .collect();
        targets.sort();
        targets.dedup();
        targets
    }
}

/// Iterator over query alignment sets with automatic backpressure.
///
/// This ensures that FastGA doesn't run ahead and fill up memory.
/// When you're processing slowly, FastGA automatically pauses.
pub struct QueryAlignmentIterator {
    receiver: Receiver<Result<QueryAlignmentSet>>,
    _producer: thread::JoinHandle<Result<()>>,
}

impl QueryAlignmentIterator {
    /// Create a new iterator from FASTA/GDB files
    ///
    /// # Arguments
    /// * `queries` - Path to query sequences (FASTA or GDB)
    /// * `targets` - Path to target database (FASTA or GDB)
    /// * `config` - Alignment configuration
    /// * `buffer_size` - How many query results to buffer (1 = strict backpressure)
    pub fn new(
        queries: impl AsRef<Path>,
        targets: impl AsRef<Path>,
        config: Config,
        buffer_size: usize,
    ) -> Result<Self> {
        let queries = queries.as_ref().to_path_buf();
        let targets = targets.as_ref().to_path_buf();

        // Create bounded channel for backpressure
        let (tx, rx) = sync_channel(buffer_size);

        // Start FastGA in background thread
        let producer = thread::spawn(move || -> Result<()> {
            // Prepare GDB files if needed
            let query_gdb = crate::ffi::prepare_database(&queries)?;
            let target_gdb = crate::ffi::prepare_database(&targets)?;

            // Run FastGA
            let mut child = Command::new(crate::embedded::get_binary_path("FastGA")?)
                .arg("-T")
                .arg(config.num_threads.to_string())
                .arg("-pafx") // PAF with extended CIGAR
                .arg("-l")
                .arg(config.min_alignment_length.to_string())
                .arg("-i")
                .arg(format!("{:.2}", config.min_identity.unwrap_or(0.7)))
                .arg(&query_gdb)
                .arg(&target_gdb)
                .stdout(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()?;

            let stdout = child.stdout.take().expect("Failed to capture stdout");
            let reader = BufReader::new(stdout);

            let mut current_query: Option<String> = None;
            let mut current_alignments: Vec<Alignment> = Vec::new();

            for line in reader.lines() {
                let line = line?;
                if line.is_empty() {
                    continue;
                }

                // Parse alignment from PAF line
                let alignment = match parse_paf_line(&line) {
                    Ok(a) => a,
                    Err(e) => {
                        eprintln!("Warning: Failed to parse PAF line: {e}");
                        continue;
                    }
                };

                // Check if we've moved to a new query
                if let Some(ref prev_query) = current_query {
                    if alignment.query_name != *prev_query {
                        // Previous query is complete - send it
                        let set = QueryAlignmentSet {
                            query_name: prev_query.clone(),
                            alignments: std::mem::take(&mut current_alignments),
                        };

                        // This blocks if receiver is busy (backpressure!)
                        if tx.send(Ok(set)).is_err() {
                            // Receiver dropped, stop
                            child.kill().ok();
                            return Ok(());
                        }
                    }
                }

                current_query = Some(alignment.query_name.clone());
                current_alignments.push(alignment);
            }

            // Send last query
            if let Some(query) = current_query {
                let set = QueryAlignmentSet {
                    query_name: query,
                    alignments: current_alignments,
                };
                let _ = tx.send(Ok(set));
            }

            child.wait()?;
            Ok(())
        });

        Ok(Self {
            receiver: rx,
            _producer: producer,
        })
    }
}

impl Iterator for QueryAlignmentIterator {
    type Item = Result<QueryAlignmentSet>;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.recv().ok()
    }
}

/// High-level API for processing queries one at a time
pub fn align_queries<F>(
    queries: impl AsRef<Path>,
    targets: impl AsRef<Path>,
    config: Config,
    mut processor: F,
) -> Result<()>
where
    F: FnMut(QueryAlignmentSet) -> Result<bool>, // Return false to stop
{
    let iter = QueryAlignmentIterator::new(queries, targets, config, 1)?;

    for result in iter {
        let set = result?;
        if !processor(set)? {
            break; // User requested stop
        }
    }

    Ok(())
}

fn parse_paf_line(line: &str) -> Result<Alignment> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return Err(anyhow::anyhow!("Invalid PAF line: too few fields"));
    }

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
        matches: fields[9].parse()?,
        block_len: fields[10].parse()?,
        mapping_quality: fields[11].parse()?,
        cigar: String::new(),
        mismatches: 0, // Will be filled from tags if present
        gap_opens: 0,  // Will be filled from tags if present
        gap_len: 0,    // Will be filled from tags if present
        tags: Vec::new(),
    };

    // Parse optional tags
    for tag in &fields[12..] {
        if let Some(cigar) = tag.strip_prefix("cg:Z:") {
            alignment.cigar = cigar.to_string();
        }
    }

    Ok(alignment)
}
