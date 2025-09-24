/// Query-by-query processor with backpressure and streaming
///
/// This allows you to process each query's alignments completely
/// before FastGA generates the next query's alignments

use std::process::{Command, Stdio, Child, ChildStdout};
use std::io::{BufReader, BufRead, Write};
use std::sync::mpsc::{sync_channel, SyncSender, Receiver};
use std::thread;
use std::path::{Path, PathBuf};
use anyhow::Result;
use tempfile::NamedTempFile;

/// Process queries one at a time with natural backpressure
pub struct QueryProcessor {
    query_gdb: PathBuf,
    target_gdb: PathBuf,
    config: crate::Config,
}

impl QueryProcessor {
    pub fn new(queries: &Path, targets: &Path, config: crate::Config) -> Result<Self> {
        Ok(Self {
            query_gdb: queries.to_path_buf(),
            target_gdb: targets.to_path_buf(),
            config,
        })
    }

    /// Method 1: Stream with backpressure using bounded channel
    /// When downstream is slow, FastGA naturally pauses
    pub fn process_with_backpressure<F>(&self, mut processor: F) -> Result<()>
    where
        F: FnMut(QueryResult) -> Result<()>,
    {
        // Bounded channel - only holds 1 query result at a time
        // This creates natural backpressure
        let (tx, rx): (SyncSender<QueryResult>, Receiver<QueryResult>) = sync_channel(1);

        // Start FastGA in a thread
        let query_gdb = self.query_gdb.clone();
        let target_gdb = self.target_gdb.clone();
        let threads = self.config.num_threads;

        let producer = thread::spawn(move || -> Result<()> {
            let mut child = Command::new("FastGA")
                .arg("-T").arg(threads.to_string())
                .arg("-pafx")
                .arg(&query_gdb)
                .arg(&target_gdb)
                .stdout(Stdio::piped())
                .spawn()?;

            let stdout = child.stdout.take().expect("Failed to capture stdout");
            let reader = BufReader::new(stdout);

            let mut current_query: Option<String> = None;
            let mut current_alignments = Vec::new();

            for line in reader.lines() {
                let line = line?;
                if line.is_empty() { continue; }

                let query_name = line.split('\t').next().unwrap_or("").to_string();

                // Detect query boundary
                if let Some(ref prev_query) = current_query {
                    if query_name != *prev_query {
                        // Send previous query result - this BLOCKS if receiver is busy!
                        let result = QueryResult {
                            query_name: prev_query.clone(),
                            alignments: std::mem::take(&mut current_alignments),
                        };

                        // This blocks until downstream consumes the result
                        if tx.send(result).is_err() {
                            // Receiver dropped, stop producing
                            child.kill()?;
                            return Ok(());
                        }
                    }
                }

                current_query = Some(query_name);
                current_alignments.push(line);
            }

            // Send last query
            if let Some(query) = current_query {
                let _ = tx.send(QueryResult {
                    query_name: query,
                    alignments: current_alignments,
                });
            }

            child.wait()?;
            Ok(())
        });

        // Process results - FastGA is paused when we're busy!
        while let Ok(result) = rx.recv() {
            processor(result)?;
        }

        producer.join().unwrap()?;
        Ok(())
    }

    /// Method 2: Write each query to a temp file
    /// Most robust for large alignments or slow downstream processing
    pub fn process_to_files<F>(&self, mut processor: F) -> Result<()>
    where
        F: FnMut(PathBuf) -> Result<()>,  // Process file path
    {
        let mut child = Command::new("FastGA")
            .arg("-T").arg(self.config.num_threads.to_string())
            .arg("-pafx")
            .arg(&self.query_gdb)
            .arg(&self.target_gdb)
            .stdout(Stdio::piped())
            .spawn()?;

        let stdout = child.stdout.take().expect("Failed to capture stdout");
        let reader = BufReader::new(stdout);

        let mut current_query: Option<String> = None;
        let mut current_file: Option<NamedTempFile> = None;

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() { continue; }

            let query_name = line.split('\t').next().unwrap_or("").to_string();

            // Detect query boundary
            if let Some(ref prev_query) = current_query {
                if query_name != *prev_query {
                    // Close current file and process it
                    if let Some(mut file) = current_file.take() {
                        file.flush()?;
                        let path = file.path().to_path_buf();

                        // Keep file alive during processing
                        let _keep = file.persist(&path)?;

                        // Process this query's file
                        processor(path.clone())?;

                        // Clean up after processing
                        std::fs::remove_file(path).ok();
                    }

                    // Start new file
                    current_file = Some(NamedTempFile::new()?);
                }
            } else {
                // First query
                current_file = Some(NamedTempFile::new()?);
            }

            // Write line to current file
            if let Some(ref mut file) = current_file {
                writeln!(file, "{}", line)?;
            }

            current_query = Some(query_name);
        }

        // Process last file
        if let Some(mut file) = current_file {
            file.flush()?;
            let path = file.path().to_path_buf();
            let _keep = file.persist(&path)?;
            processor(path.clone())?;
            std::fs::remove_file(path).ok();
        }

        child.wait()?;
        Ok(())
    }

    /// Method 3: Direct pipe to downstream command
    /// Most efficient for CLI tools that accept PAF on stdin
    pub fn pipe_to_command(&self, downstream_cmd: &str, downstream_args: &[&str]) -> Result<()> {
        let mut fastga = Command::new("FastGA")
            .arg("-T").arg(self.config.num_threads.to_string())
            .arg("-pafx")
            .arg(&self.query_gdb)
            .arg(&self.target_gdb)
            .stdout(Stdio::piped())
            .spawn()?;

        let mut downstream = Command::new(downstream_cmd)
            .args(downstream_args)
            .stdin(Stdio::piped())
            .spawn()?;

        // Connect FastGA stdout to downstream stdin
        if let Some(fastga_out) = fastga.stdout.take() {
            if let Some(mut downstream_in) = downstream.stdin.take() {
                // This creates a direct pipe - very efficient!
                std::io::copy(&mut &fastga_out, &mut downstream_in)?;
            }
        }

        fastga.wait()?;
        downstream.wait()?;
        Ok(())
    }

    /// Method 4: Process with ability to pause/resume
    /// Uses Unix signals to pause FastGA when needed
    #[cfg(unix)]
    pub fn process_with_pause<F>(&self, mut processor: F) -> Result<()>
    where
        F: FnMut(QueryResult, &mut PauseHandle) -> Result<()>,
    {
        use nix::sys::signal::{kill, Signal};
        use nix::unistd::Pid;

        let mut child = Command::new("FastGA")
            .arg("-T").arg(self.config.num_threads.to_string())
            .arg("-pafx")
            .arg(&self.query_gdb)
            .arg(&self.target_gdb)
            .stdout(Stdio::piped())
            .spawn()?;

        let pid = Pid::from_raw(child.id() as i32);
        let mut pause_handle = PauseHandle { pid };

        let stdout = child.stdout.take().expect("Failed to capture stdout");
        let reader = BufReader::new(stdout);

        let mut current_query: Option<String> = None;
        let mut current_alignments = Vec::new();

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() { continue; }

            let query_name = line.split('\t').next().unwrap_or("").to_string();

            // Detect query boundary
            if let Some(ref prev_query) = current_query {
                if query_name != *prev_query {
                    // Process complete query
                    let result = QueryResult {
                        query_name: prev_query.clone(),
                        alignments: std::mem::take(&mut current_alignments),
                    };

                    // Processor can pause/resume FastGA as needed
                    processor(result, &mut pause_handle)?;
                }
            }

            current_query = Some(query_name);
            current_alignments.push(line);
        }

        // Process last query
        if let Some(query) = current_query {
            processor(QueryResult {
                query_name: query,
                alignments: current_alignments,
            }, &mut pause_handle)?;
        }

        child.wait()?;
        Ok(())
    }
}

pub struct QueryResult {
    pub query_name: String,
    pub alignments: Vec<String>,  // PAF lines
}

#[cfg(unix)]
pub struct PauseHandle {
    pid: nix::unistd::Pid,
}

#[cfg(unix)]
impl PauseHandle {
    pub fn pause(&self) -> Result<()> {
        use nix::sys::signal::{kill, Signal};
        kill(self.pid, Signal::SIGSTOP)?;
        Ok(())
    }

    pub fn resume(&self) -> Result<()> {
        use nix::sys::signal::{kill, Signal};
        kill(self.pid, Signal::SIGCONT)?;
        Ok(())
    }
}

/// Example usage patterns:
///
/// ```rust
/// // Pattern 1: Natural backpressure
/// processor.process_with_backpressure(|result| {
///     // FastGA is automatically paused while we process this query
///     expensive_computation(&result)?;
///     Ok(())
/// })?;
///
/// // Pattern 2: Write to files for another tool
/// processor.process_to_files(|paf_file| {
///     Command::new("downstream_tool")
///         .arg(paf_file)
///         .status()?;
///     Ok(())
/// })?;
///
/// // Pattern 3: Direct pipe (most efficient)
/// processor.pipe_to_command("sort", &["-k1,1", "-k3,3n"])?;
///
/// // Pattern 4: Explicit pause control (Unix only)
/// processor.process_with_pause(|result, pause| {
///     if result.alignments.len() > 1000 {
///         pause.pause()?;  // Stop FastGA
///         process_large_result(&result)?;
///         pause.resume()?; // Continue FastGA
///     }
///     Ok(())
/// })?;
/// ```