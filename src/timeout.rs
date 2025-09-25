//! Timeout support for FastGA alignment operations

use crate::error::{Result, FastGAError};
use crate::{Alignment, Alignments, Config, FastGA};
use std::path::Path;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

/// Extension trait to add timeout support to FastGA
pub trait TimeoutExt {
    /// Align files with a timeout
    fn align_files_timeout(
        &self,
        genome1: &Path,
        genome2: &Path,
        timeout: Duration,
    ) -> Result<Alignments>;
}

impl TimeoutExt for FastGA {
    fn align_files_timeout(
        &self,
        genome1: &Path,
        genome2: &Path,
        timeout: Duration,
    ) -> Result<Alignments> {
        let genome1 = genome1.to_path_buf();
        let genome2 = genome2.to_path_buf();
        let inner_clone = self.inner.clone();

        let (tx, rx) = mpsc::channel();

        thread::spawn(move || {
            let result = inner_clone.align_files(&genome1, &genome2);
            let _ = tx.send(result);
        });

        match rx.recv_timeout(timeout) {
            Ok(result) => result,
            Err(mpsc::RecvTimeoutError::Timeout) => {
                Err(FastGAError::Other("Alignment operation timed out".to_string()))
            }
            Err(mpsc::RecvTimeoutError::Disconnected) => {
                Err(FastGAError::Other("Alignment thread crashed".to_string()))
            }
        }
    }
}

/// Wrapper with built-in timeout and progress support
pub struct TimeoutAligner {
    config: Config,
    timeout: Option<Duration>,
    progress_callback: Option<Box<dyn Fn(&str, &str) + Send + Sync>>,
}

impl TimeoutAligner {
    pub fn new(config: Config) -> Self {
        Self {
            config,
            timeout: None,
            progress_callback: None,
        }
    }

    pub fn with_timeout(mut self, timeout: Duration) -> Self {
        self.timeout = Some(timeout);
        self
    }

    pub fn with_progress<F>(mut self, callback: F) -> Self
    where
        F: Fn(&str, &str) + Send + Sync + 'static,
    {
        self.progress_callback = Some(Box::new(callback));
        self
    }

    pub fn align_files(&self, genome1: &Path, genome2: &Path) -> Result<Alignments> {
        if let Some(ref callback) = self.progress_callback {
            callback("starting", "Initializing alignment");
        }

        let aligner = FastGA::new(self.config.clone())?;

        if let Some(ref callback) = self.progress_callback {
            callback("aligning", "Running FastGA alignment");
        }

        let result = if let Some(timeout) = self.timeout {
            aligner.align_files_timeout(genome1, genome2, timeout)
        } else {
            aligner.align_files(genome1, genome2)
        };

        if let Some(ref callback) = self.progress_callback {
            match &result {
                Ok(alignments) => {
                    callback(
                        "complete",
                        &format!("Alignment complete: {} alignments found", alignments.len()),
                    );
                }
                Err(e) => {
                    callback("error", &format!("Alignment failed: {}", e));
                }
            }
        }

        result
    }
}