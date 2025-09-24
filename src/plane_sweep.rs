//! Plane sweep filtering for streaming alignments
//!
//! This module provides immediate plane sweep filtering per query to handle
//! repetitive genomes efficiently without writing massive intermediate files.

use crate::{Alignment, Result};
use std::collections::{BTreeMap, HashMap};
use std::cmp::Ordering;

/// Configuration for plane sweep filtering
#[derive(Debug, Clone)]
pub struct PlaneSweepConfig {
    /// Maximum alignments to keep per query (0 = unlimited)
    pub max_per_query: usize,
    /// Maximum alignments to keep per query-target pair
    pub max_per_target: usize,
    /// Minimum identity threshold
    pub min_identity: f64,
    /// Minimum alignment length
    pub min_length: usize,
    /// Overlap fraction threshold (0.0-1.0)
    pub max_overlap: f64,
}

impl Default for PlaneSweepConfig {
    fn default() -> Self {
        PlaneSweepConfig {
            max_per_query: 100,      // Keep top 100 per query
            max_per_target: 1,       // 1:1 mapping by default
            min_identity: 0.0,
            min_length: 0,
            max_overlap: 0.5,        // Filter if >50% overlap with better alignment
        }
    }
}

/// Plane sweep filter that processes alignments per query
pub struct PlaneSweepFilter {
    config: PlaneSweepConfig,
    current_query: Option<String>,
    active_alignments: BTreeMap<OrderedAlignment, Alignment>,
    kept_alignments: Vec<Alignment>,
    stats: FilterStats,
}

/// Statistics from filtering
#[derive(Debug, Default, Clone)]
pub struct FilterStats {
    pub total_processed: usize,
    pub kept: usize,
    pub filtered_by_identity: usize,
    pub filtered_by_length: usize,
    pub filtered_by_overlap: usize,
    pub filtered_by_limit: usize,
}

/// Wrapper for alignment ordering by score
#[derive(Debug, Clone)]
struct OrderedAlignment {
    score: OrderedFloat,
    rank: usize,
}

#[derive(Debug, Clone, PartialEq)]
struct OrderedFloat(f64);

impl Eq for OrderedFloat {}

impl PartialOrd for OrderedFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Higher scores first
        other.0.partial_cmp(&self.0)
    }
}

impl Ord for OrderedFloat {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

impl PartialEq for OrderedAlignment {
    fn eq(&self, other: &Self) -> bool {
        self.rank == other.rank
    }
}

impl Eq for OrderedAlignment {}

impl PartialOrd for OrderedAlignment {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedAlignment {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.score.cmp(&other.score) {
            Ordering::Equal => self.rank.cmp(&other.rank),
            other => other,
        }
    }
}

impl PlaneSweepFilter {
    /// Create a new plane sweep filter
    pub fn new(config: PlaneSweepConfig) -> Self {
        PlaneSweepFilter {
            config,
            current_query: None,
            active_alignments: BTreeMap::new(),
            kept_alignments: Vec::new(),
            stats: FilterStats::default(),
        }
    }

    /// Process a single alignment
    /// Returns true if alignment should be kept
    pub fn process(&mut self, alignment: Alignment) -> bool {
        self.stats.total_processed += 1;

        // Check if we've moved to a new query
        if self.current_query.as_ref() != Some(&alignment.query_name) {
            self.flush_query();
            self.current_query = Some(alignment.query_name.clone());
        }

        // Apply basic filters
        if alignment.identity() < self.config.min_identity {
            self.stats.filtered_by_identity += 1;
            return false;
        }

        let length = alignment.query_end - alignment.query_start;
        if length < self.config.min_length {
            self.stats.filtered_by_length += 1;
            return false;
        }

        // Calculate plane sweep score
        let score = alignment.identity() * (length as f64).ln();

        // Check overlap with existing alignments
        if self.has_significant_overlap(&alignment) {
            self.stats.filtered_by_overlap += 1;
            return false;
        }

        // Add to active set
        let ordered = OrderedAlignment {
            score: OrderedFloat(score),
            rank: self.stats.total_processed,
        };

        self.active_alignments.insert(ordered, alignment);

        // Enforce per-query limit
        if self.config.max_per_query > 0 {
            while self.active_alignments.len() > self.config.max_per_query {
                self.active_alignments.pop_last();
                self.stats.filtered_by_limit += 1;
            }
        }

        true
    }

    /// Check if alignment significantly overlaps with better alignments
    fn has_significant_overlap(&self, alignment: &Alignment) -> bool {
        for (_, existing) in self.active_alignments.iter() {
            if existing.target_name != alignment.target_name {
                continue;
            }

            let overlap = Self::calculate_overlap(
                existing.query_start, existing.query_end,
                alignment.query_start, alignment.query_end,
            );

            let alignment_length = alignment.query_end - alignment.query_start;
            let overlap_fraction = overlap as f64 / alignment_length as f64;

            if overlap_fraction > self.config.max_overlap {
                return true;
            }
        }
        false
    }

    /// Calculate overlap between two intervals
    fn calculate_overlap(start1: usize, end1: usize, start2: usize, end2: usize) -> usize {
        let overlap_start = start1.max(start2);
        let overlap_end = end1.min(end2);
        if overlap_end > overlap_start {
            overlap_end - overlap_start
        } else {
            0
        }
    }

    /// Flush alignments for current query
    fn flush_query(&mut self) {
        // Move active alignments to kept list
        let alignments = std::mem::take(&mut self.active_alignments);
        for (_, alignment) in alignments {
            self.kept_alignments.push(alignment);
            self.stats.kept += 1;
        }
    }

    /// Finish processing and return kept alignments
    pub fn finish(mut self) -> (Vec<Alignment>, FilterStats) {
        self.flush_query();
        (self.kept_alignments, self.stats)
    }

    /// Get current statistics
    pub fn stats(&self) -> &FilterStats {
        &self.stats
    }
}

/// Per-query streaming processor that ensures fair mapping distribution
pub struct QueryWiseSweep {
    config: PlaneSweepConfig,
    per_target_counts: HashMap<String, HashMap<String, usize>>,
}

impl QueryWiseSweep {
    /// Create a new query-wise sweep processor
    pub fn new(config: PlaneSweepConfig) -> Self {
        QueryWiseSweep {
            config,
            per_target_counts: HashMap::new(),
        }
    }

    /// Process alignments for a single query against all targets
    /// This ensures each query gets fair representation
    pub fn process_query<I>(&mut self, query_name: &str, alignments: I) -> Vec<Alignment>
    where
        I: Iterator<Item = Alignment>,
    {
        let mut filter = PlaneSweepFilter::new(self.config.clone());
        let mut query_results = Vec::new();

        // Group by target to enforce per-target limits
        let mut by_target: HashMap<String, Vec<Alignment>> = HashMap::new();

        for alignment in alignments {
            if filter.process(alignment.clone()) {
                by_target
                    .entry(alignment.target_name.clone())
                    .or_insert_with(Vec::new)
                    .push(alignment);
            }
        }

        // Apply per-target limits
        for (target, mut target_alignments) in by_target {
            // Sort by score
            target_alignments.sort_by(|a, b| {
                let score_a = a.identity() * ((a.query_end - a.query_start) as f64).ln();
                let score_b = b.identity() * ((b.query_end - b.query_start) as f64).ln();
                score_b.partial_cmp(&score_a).unwrap_or(Ordering::Equal)
            });

            // Keep only top N per target
            let limit = self.config.max_per_target;
            for alignment in target_alignments.into_iter().take(limit) {
                query_results.push(alignment);

                // Track per-target counts for statistics
                *self.per_target_counts
                    .entry(query_name.to_string())
                    .or_insert_with(HashMap::new)
                    .entry(target.clone())
                    .or_insert(0) += 1;
            }
        }

        query_results
    }

    /// Get distribution statistics
    pub fn get_distribution_stats(&self) -> String {
        let mut stats = String::new();
        for (query, targets) in &self.per_target_counts {
            let total = targets.values().sum::<usize>();
            let unique_targets = targets.len();
            stats.push_str(&format!(
                "{}: {} alignments to {} unique targets\n",
                query, total, unique_targets
            ));
        }
        stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plane_sweep_scoring() {
        let config = PlaneSweepConfig::default();
        let mut filter = PlaneSweepFilter::new(config);

        let alignment = Alignment {
            query_name: "query1".to_string(),
            query_len: 1000,
            query_start: 0,
            query_end: 100,
            strand: '+',
            target_name: "target1".to_string(),
            target_len: 1000,
            target_start: 0,
            target_end: 100,
            matches: 95,
            block_len: 100,
            mapping_quality: 255,
            cigar: "95=5X".to_string(),
            mismatches: 5,
            gap_opens: 0,
            gap_len: 0,
            tags: Vec::new(),
        };

        assert!(filter.process(alignment));

        let (kept, stats) = filter.finish();
        assert_eq!(kept.len(), 1);
        assert_eq!(stats.kept, 1);
    }

    #[test]
    fn test_overlap_filtering() {
        let mut config = PlaneSweepConfig::default();
        config.max_overlap = 0.5;
        let mut filter = PlaneSweepFilter::new(config);

        // First alignment
        let alignment1 = Alignment {
            query_name: "query1".to_string(),
            query_len: 1000,
            query_start: 0,
            query_end: 100,
            strand: '+',
            target_name: "target1".to_string(),
            target_len: 1000,
            target_start: 0,
            target_end: 100,
            matches: 100,
            block_len: 100,
            mapping_quality: 255,
            cigar: "100=".to_string(),
            mismatches: 0,
            gap_opens: 0,
            gap_len: 0,
            tags: Vec::new(),
        };

        // Overlapping alignment (60% overlap)
        let alignment2 = Alignment {
            query_name: "query1".to_string(),
            query_len: 1000,
            query_start: 40,
            query_end: 140,
            strand: '+',
            target_name: "target1".to_string(),
            target_len: 1000,
            target_start: 40,
            target_end: 140,
            matches: 90,
            block_len: 100,
            mapping_quality: 255,
            cigar: "90=10X".to_string(),
            mismatches: 10,
            gap_opens: 0,
            gap_len: 0,
            tags: Vec::new(),
        };

        assert!(filter.process(alignment1));
        assert!(!filter.process(alignment2)); // Should be filtered due to overlap

        let (kept, stats) = filter.finish();
        assert_eq!(kept.len(), 1);
        assert_eq!(stats.filtered_by_overlap, 1);
    }
}