//! Integration module for embedding FastGA-rs into SweepGA
//!
//! This module provides the API for SweepGA to directly call FastGA
//! and process alignments without intermediate files.

use crate::{Config, Alignment, Result};
use crate::plane_sweep::PlaneSweepConfig;
use std::path::Path;

/// Direct integration point for SweepGA
///
/// This struct is designed to be embedded into SweepGA's filtering pipeline,
/// allowing it to generate alignments on-demand and apply plane sweep filtering
/// immediately, avoiding massive intermediate files for repetitive genomes.
pub struct IntegratedAligner {
    align_config: Config,
    sweep_config: Option<PlaneSweepConfig>,
}

impl IntegratedAligner {
    /// Create a new integrated aligner for use within SweepGA
    pub fn new(align_config: Config) -> Self {
        IntegratedAligner {
            align_config,
            sweep_config: None,
        }
    }

    /// Enable plane sweep pre-filtering to handle repetitive genomes
    pub fn with_plane_sweep(mut self, config: PlaneSweepConfig) -> Self {
        self.sweep_config = Some(config);
        self
    }

    /// Process a single query contig against all target contigs
    ///
    /// This is the key method for SweepGA integration:
    /// - Takes one query contig at a time
    /// - Aligns against all targets
    /// - Applies immediate plane sweep filtering if configured
    /// - Returns only the essential alignments
    ///
    /// This allows SweepGA to process massive genomes without disk explosion
    pub fn align_query_contig<F>(
        &self,
        query_name: &str,
        query_seq: &[u8],
        target_file: &Path,
        mut callback: F,
    ) -> Result<Vec<Alignment>>
    where
        F: FnMut(&Alignment) -> bool,
    {
        // Create temp file for single query
        let temp_dir = tempfile::TempDir::new()?;
        let query_path = temp_dir.path().join("query.fa");
        std::fs::write(&query_path, format!(
            ">{}\n{}\n",
            query_name,
            std::str::from_utf8(query_seq).unwrap_or("")
        ))?;

        // Run FastGA on this single query
        let aligner = crate::FastGA::new(self.align_config.clone())?;
        let raw_alignments = aligner.align_files(&query_path, target_file)?;

        // Apply filtering
        let mut kept = Vec::new();
        for alignment in raw_alignments.alignments {
            // Apply plane sweep filter if configured
            if let Some(ref sweep_config) = self.sweep_config {
                if alignment.identity() < sweep_config.min_identity {
                    continue;
                }
                if alignment.query_end - alignment.query_start < sweep_config.min_length {
                    continue;
                }
                // Check if we've hit the per-query limit
                if sweep_config.max_per_query > 0 && kept.len() >= sweep_config.max_per_query {
                    break;
                }
            }

            // Let SweepGA's callback decide too
            if callback(&alignment) {
                kept.push(alignment);
            }
        }

        Ok(kept)
    }

    /// Process all queries from a FASTA file against targets
    ///
    /// This method is designed for SweepGA's batch processing:
    /// - Reads queries one by one from a FASTA file
    /// - Processes each independently (fair distribution)
    /// - Applies both plane sweep and SweepGA's filters
    /// - Minimizes memory usage for repetitive genomes
    pub fn process_all_queries<F, G>(
        &self,
        query_file: &Path,
        target_file: &Path,
        mut per_alignment_callback: F,
        mut per_query_callback: G,
    ) -> Result<usize>
    where
        F: FnMut(&Alignment) -> bool,
        G: FnMut(&str, Vec<Alignment>),
    {
        use std::io::{BufRead, BufReader};

        let file = std::fs::File::open(query_file)?;
        let reader = BufReader::new(file);

        let mut current_name = String::new();
        let mut current_seq = Vec::new();
        let mut total_kept = 0;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                // Process previous sequence if exists
                if !current_name.is_empty() {
                    let alignments = self.align_query_contig(
                        &current_name,
                        &current_seq,
                        target_file,
                        &mut per_alignment_callback,
                    )?;

                    total_kept += alignments.len();
                    per_query_callback(&current_name, alignments);
                }

                // Start new sequence
                current_name = line[1..].split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                current_seq.clear();
            } else {
                current_seq.extend_from_slice(line.as_bytes());
            }
        }

        // Process last sequence
        if !current_name.is_empty() {
            let alignments = self.align_query_contig(
                &current_name,
                &current_seq,
                target_file,
                &mut per_alignment_callback,
            )?;

            total_kept += alignments.len();
            per_query_callback(&current_name, alignments);
        }

        Ok(total_kept)
    }
}

/// Streaming iterator for memory-efficient processing
///
/// This is designed for SweepGA to process alignments as a stream
/// without holding everything in memory at once.
pub struct AlignmentStream {
    aligner: IntegratedAligner,
    queries: Vec<(String, Vec<u8>)>,
    target_file: std::path::PathBuf,
    current_query_idx: usize,
    current_alignments: std::vec::IntoIter<Alignment>,
}

impl AlignmentStream {
    /// Create a new alignment stream
    pub fn new(
        aligner: IntegratedAligner,
        query_file: &Path,
        target_file: &Path,
    ) -> Result<Self> {
        // Load queries (in practice, could stream these too)
        let queries = load_fasta_sequences(query_file)?;

        Ok(AlignmentStream {
            aligner,
            queries,
            target_file: target_file.to_path_buf(),
            current_query_idx: 0,
            current_alignments: Vec::new().into_iter(),
        })
    }
}

impl Iterator for AlignmentStream {
    type Item = Alignment;

    fn next(&mut self) -> Option<Self::Item> {
        // Try to get next alignment from current query
        if let Some(alignment) = self.current_alignments.next() {
            return Some(alignment);
        }

        // Move to next query
        while self.current_query_idx < self.queries.len() {
            let (query_name, query_seq) = &self.queries[self.current_query_idx];
            self.current_query_idx += 1;

            // Align this query
            if let Ok(alignments) = self.aligner.align_query_contig(
                query_name,
                query_seq,
                &self.target_file,
                |_| true, // Accept all for now
            ) {
                self.current_alignments = alignments.into_iter();
                if let Some(alignment) = self.current_alignments.next() {
                    return Some(alignment);
                }
            }
        }

        None
    }
}

// Helper function to load FASTA sequences
fn load_fasta_sequences(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_name.is_empty() {
                sequences.push((current_name.clone(), current_seq.clone()));
            }
            current_name = line[1..].split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_seq.clear();
        } else {
            current_seq.extend_from_slice(line.as_bytes());
        }
    }

    if !current_name.is_empty() {
        sequences.push((current_name, current_seq));
    }

    Ok(sequences)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrated_aligner_creation() {
        let aligner = IntegratedAligner::new(Config::default())
            .with_plane_sweep(PlaneSweepConfig {
                max_per_query: 100,
                max_per_target: 1,
                min_identity: 0.9,
                min_length: 1000,
                max_overlap: 0.5,
            });

        assert!(aligner.sweep_config.is_some());
    }
}