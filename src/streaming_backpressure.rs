/// Clean streaming API with natural backpressure
/// No temp files - just memory-bounded streaming

use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead};
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread;
use std::path::Path;
use anyhow::Result;
use crate::{Alignment, Config};

/// Stream alignments with automatic backpressure
/// When you're slow processing a query, FastGA naturally pauses
pub struct BackpressureStream {
    receiver: Receiver<QueryBatch>,
    _producer: thread::JoinHandle<Result<()>>,
}

impl BackpressureStream {
    /// Start FastGA and return a stream of query results
    /// Buffer size controls how many queries can be queued
    pub fn new(
        query_path: &Path,
        target_path: &Path,
        config: Config,
        buffer_size: usize,  // How many queries to buffer (1 = strict backpressure)
    ) -> Result<Self> {
        // Convert to GDB if needed
        let query_gdb = crate::ffi::prepare_database(query_path)?;
        let target_gdb = crate::ffi::prepare_database(target_path)?;

        // Bounded channel creates backpressure
        let (tx, rx) = sync_channel(buffer_size);

        // Start FastGA in background thread
        let producer = thread::spawn(move || -> Result<()> {
            let mut child = Command::new("FastGA")
                .arg("-T").arg(config.num_threads.to_string())
                .arg("-pafx")  // PAF with extended CIGAR
                .arg(&query_gdb)
                .arg(&target_gdb)
                .stdout(Stdio::piped())
                .stderr(Stdio::null())  // Suppress progress output
                .spawn()?;

            let stdout = child.stdout.take().expect("Failed to capture stdout");
            let reader = BufReader::new(stdout);

            let mut current_query: Option<String> = None;
            let mut current_alignments = Vec::new();

            for line in reader.lines() {
                let line = line?;
                if line.is_empty() { continue; }

                // Parse query name from PAF
                let query_name = line.split('\t').next()
                    .ok_or_else(|| anyhow::anyhow!("Invalid PAF line"))?
                    .to_string();

                // Check for query boundary
                if let Some(ref prev_query) = current_query {
                    if query_name != *prev_query {
                        // Send complete query batch
                        // This BLOCKS if receiver is busy (backpressure!)
                        let batch = QueryBatch {
                            query_name: prev_query.clone(),
                            alignments: parse_alignments(&current_alignments)?,
                        };

                        if tx.send(batch).is_err() {
                            // Receiver dropped, stop
                            child.kill().ok();
                            return Ok(());
                        }

                        current_alignments.clear();
                    }
                }

                current_query = Some(query_name);
                current_alignments.push(line);
            }

            // Send last query
            if let Some(query) = current_query {
                let batch = QueryBatch {
                    query_name: query,
                    alignments: parse_alignments(&current_alignments)?,
                };
                let _ = tx.send(batch);
            }

            child.wait()?;
            Ok(())
        });

        Ok(Self {
            receiver: rx,
            _producer: producer,
        })
    }

    /// Get next query's alignments (blocks until available)
    pub fn next_query(&self) -> Option<QueryBatch> {
        self.receiver.recv().ok()
    }

    /// Iterate over queries
    pub fn iter(&self) -> QueryIterator {
        QueryIterator { stream: self }
    }
}

pub struct QueryIterator<'a> {
    stream: &'a BackpressureStream,
}

impl<'a> Iterator for QueryIterator<'a> {
    type Item = QueryBatch;

    fn next(&mut self) -> Option<Self::Item> {
        self.stream.next_query()
    }
}

/// A complete set of alignments for one query
pub struct QueryBatch {
    pub query_name: String,
    pub alignments: Vec<Alignment>,
}

impl QueryBatch {
    /// Write this query's alignments to a file if needed
    pub fn write_paf(&self, path: &Path) -> Result<()> {
        use std::fs::File;
        use std::io::Write;

        let mut file = File::create(path)?;
        for aln in &self.alignments {
            writeln!(file, "{}", aln.to_paf_line())?;
        }
        Ok(())
    }

    /// Get total bases aligned for this query
    pub fn total_aligned_bases(&self) -> usize {
        self.alignments.iter()
            .map(|a| a.query_end - a.query_start)
            .sum()
    }

    /// Filter alignments by identity
    pub fn filter_by_identity(&mut self, min_identity: f64) {
        self.alignments.retain(|a| a.identity() >= min_identity);
    }
}

fn parse_alignments(paf_lines: &[String]) -> Result<Vec<Alignment>> {
    paf_lines.iter()
        .map(|line| parse_paf_line(line))
        .collect()
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
        mismatches: 0,
        gap_opens: 0,
        gap_len: 0,
        tags: Vec::new(),
    };

    // Parse CIGAR from tags
    for tag in &fields[12..] {
        if tag.starts_with("cg:Z:") {
            alignment.cigar = tag[5..].to_string();
        } else if tag.starts_with("NM:i:") {
            // Parse number of mismatches
            if let Ok(nm) = tag[5..].parse::<usize>() {
                alignment.mismatches = nm;
            }
        }
    }

    Ok(alignment)
}

/// Example usage:
/// ```rust
/// // Create stream with buffer for 1 query (strict backpressure)
/// let stream = BackpressureStream::new(
///     "queries.fa",
///     "targets.fa",
///     Config::default(),
///     1  // Only buffer 1 query ahead
/// )?;
///
/// // Process each query as it completes
/// for batch in stream.iter() {
///     println!("Processing {} with {} alignments",
///              batch.query_name, batch.alignments.len());
///
///     // FastGA is paused here if we're slow!
///     expensive_computation(&batch)?;
///
///     // Or write to file for another tool
///     batch.write_paf(&format!("{}.paf", batch.query_name))?;
///
///     // Or pipe to another command
///     let filtered: Vec<_> = batch.alignments.iter()
///         .filter(|a| a.identity() > 0.95)
///         .collect();
/// }
/// ```
///
/// Benefits:
/// - No temp files
/// - Natural backpressure (FastGA pauses when you're slow)
/// - Memory bounded (only buffers N queries)
/// - Clean iteration API