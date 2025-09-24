//! Alignment representation and format conversion.
//!
//! This module provides structures for representing alignments and converting
//! between different alignment formats (PAF, PSL, trace points).

use crate::error::{Result, FastGAError};
use std::fmt;

/// A single local alignment between two sequences.
///
/// This structure represents one alignment with full CIGAR information
/// using extended CIGAR format (with '=' and 'X' operators).
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Query sequence name
    pub query_name: String,

    /// Query sequence length
    pub query_len: usize,

    /// Query start position (0-based)
    pub query_start: usize,

    /// Query end position (0-based, exclusive)
    pub query_end: usize,

    /// Strand ('+' for forward, '-' for reverse)
    pub strand: char,

    /// Target sequence name
    pub target_name: String,

    /// Target sequence length
    pub target_len: usize,

    /// Target start position (0-based)
    pub target_start: usize,

    /// Target end position (0-based, exclusive)
    pub target_end: usize,

    /// Number of matching bases
    pub matches: usize,

    /// Alignment block length
    pub block_len: usize,

    /// Mapping quality (0-255)
    pub mapping_quality: u8,

    /// Extended CIGAR string with '=' and 'X' operators
    pub cigar: String,

    /// Number of mismatches
    pub mismatches: usize,

    /// Number of gap opens
    pub gap_opens: usize,

    /// Total gap length
    pub gap_len: usize,

    /// Additional SAM-style tags
    pub tags: Vec<(String, String)>,
}

impl Alignment {
    /// Calculates the identity fraction of the alignment.
    ///
    /// Identity is calculated as matches / (matches + mismatches + gaps)
    pub fn identity(&self) -> f64 {
        let total = self.matches + self.mismatches + self.gap_len;
        if total > 0 {
            self.matches as f64 / total as f64
        } else {
            0.0
        }
    }

    /// Returns the query coverage fraction.
    pub fn query_coverage(&self) -> f64 {
        (self.query_end - self.query_start) as f64 / self.query_len as f64
    }

    /// Returns the target coverage fraction.
    pub fn target_coverage(&self) -> f64 {
        (self.target_end - self.target_start) as f64 / self.target_len as f64
    }

    /// Formats the alignment as a PAF line.
    pub fn to_paf(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
            self.query_name,
            self.query_len,
            self.query_start,
            self.query_end,
            self.strand,
            self.target_name,
            self.target_len,
            self.target_start,
            self.target_end,
            self.matches,
            self.block_len,
            self.mapping_quality,
            self.cigar
        )
    }

    /// Parses a PAF line into an Alignment.
    pub fn from_paf_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            return Err(FastGAError::PafParseError(
                format!("PAF line has {} fields, expected at least 12", fields.len())
            ));
        }

        let mut alignment = Alignment {
            query_name: fields[0].to_string(),
            query_len: fields[1].parse().map_err(|_|
                FastGAError::PafParseError("Invalid query length".to_string()))?,
            query_start: fields[2].parse().map_err(|_|
                FastGAError::PafParseError("Invalid query start".to_string()))?,
            query_end: fields[3].parse().map_err(|_|
                FastGAError::PafParseError("Invalid query end".to_string()))?,
            strand: fields[4].chars().next().ok_or_else(||
                FastGAError::PafParseError("Invalid strand".to_string()))?,
            target_name: fields[5].to_string(),
            target_len: fields[6].parse().map_err(|_|
                FastGAError::PafParseError("Invalid target length".to_string()))?,
            target_start: fields[7].parse().map_err(|_|
                FastGAError::PafParseError("Invalid target start".to_string()))?,
            target_end: fields[8].parse().map_err(|_|
                FastGAError::PafParseError("Invalid target end".to_string()))?,
            matches: fields[9].parse().map_err(|_|
                FastGAError::PafParseError("Invalid matches".to_string()))?,
            block_len: fields[10].parse().map_err(|_|
                FastGAError::PafParseError("Invalid block length".to_string()))?,
            mapping_quality: fields[11].parse().map_err(|_|
                FastGAError::PafParseError("Invalid mapping quality".to_string()))?,
            cigar: String::new(),
            mismatches: 0,
            gap_opens: 0,
            gap_len: 0,
            tags: Vec::new(),
        };

        // Parse optional tags
        for field in &fields[12..] {
            if let Some((tag, value)) = field.split_once(':') {
                if tag == "cg" && value.starts_with("Z:") {
                    alignment.cigar = value[2..].to_string();
                    alignment.parse_cigar_stats()?;
                } else {
                    alignment.tags.push((tag.to_string(), value.to_string()));
                }
            }
        }

        Ok(alignment)
    }

    /// Parse CIGAR string to extract mismatch and gap statistics.
    fn parse_cigar_stats(&mut self) -> Result<()> {
        let mut mismatches = 0;
        let mut gap_opens = 0;
        let mut gap_len = 0;
        let mut num_str = String::new();
        let mut in_gap = false;

        for ch in self.cigar.chars() {
            if ch.is_ascii_digit() {
                num_str.push(ch);
            } else {
                let count: usize = num_str.parse().map_err(|_|
                    FastGAError::CigarParseError("Invalid number in CIGAR".to_string()))?;
                num_str.clear();

                match ch {
                    'X' => mismatches += count,  // Mismatch
                    'I' | 'D' => {  // Insertion or Deletion
                        gap_len += count;
                        if !in_gap {
                            gap_opens += 1;
                            in_gap = true;
                        }
                    }
                    '=' | 'M' => in_gap = false,  // Match
                    _ => {}
                }
            }
        }

        self.mismatches = mismatches;
        self.gap_opens = gap_opens;
        self.gap_len = gap_len;

        Ok(())
    }
}

/// Collection of alignments with format conversion capabilities.
#[derive(Debug)]
pub struct Alignments {
    /// Vector of individual alignments
    pub alignments: Vec<Alignment>,
}

impl Alignments {
    /// Creates a new empty alignment collection.
    pub fn new() -> Self {
        Alignments {
            alignments: Vec::new(),
        }
    }

    /// Creates alignments from PAF format text.
    pub fn from_paf(paf_text: &str) -> Result<Self> {
        let mut alignments = Vec::new();

        for line in paf_text.lines() {
            if !line.is_empty() && !line.starts_with('#') {
                alignments.push(Alignment::from_paf_line(line)?);
            }
        }

        Ok(Alignments { alignments })
    }

    /// Converts all alignments to PAF format.
    pub fn to_paf(&self) -> Result<String> {
        let lines: Vec<String> = self.alignments
            .iter()
            .map(|a| a.to_paf())
            .collect();
        Ok(lines.join("\n"))
    }

    /// Returns an iterator over the alignments.
    pub fn iter(&self) -> impl Iterator<Item = &Alignment> {
        self.alignments.iter()
    }

    /// Returns the number of alignments.
    pub fn len(&self) -> usize {
        self.alignments.len()
    }

    /// Returns true if there are no alignments.
    pub fn is_empty(&self) -> bool {
        self.alignments.is_empty()
    }

    /// Filters alignments by a predicate.
    pub fn filter<F>(&mut self, predicate: F)
    where
        F: Fn(&Alignment) -> bool,
    {
        self.alignments.retain(predicate);
    }

    /// Sorts alignments by query position.
    pub fn sort_by_query(&mut self) {
        self.alignments.sort_by_key(|a| (a.query_name.clone(), a.query_start));
    }

    /// Sorts alignments by target position.
    pub fn sort_by_target(&mut self) {
        self.alignments.sort_by_key(|a| (a.target_name.clone(), a.target_start));
    }

    /// Groups alignments by query sequence.
    pub fn group_by_query(&self) -> std::collections::HashMap<String, Vec<&Alignment>> {
        let mut groups = std::collections::HashMap::new();

        for alignment in &self.alignments {
            groups.entry(alignment.query_name.clone())
                .or_insert_with(Vec::new)
                .push(alignment);
        }

        groups
    }
}

impl fmt::Display for Alignments {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Alignments: {} total", self.alignments.len())?;
        for (i, alignment) in self.alignments.iter().enumerate() {
            writeln!(
                f,
                "[{}] {} -> {}: {:.1}% identity, {} bp",
                i,
                alignment.query_name,
                alignment.target_name,
                alignment.identity() * 100.0,
                alignment.query_end - alignment.query_start
            )?;
        }
        Ok(())
    }
}