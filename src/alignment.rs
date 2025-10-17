//! Alignment representation and format conversion.
//!
//! This module provides structures for representing alignments and converting
//! between different alignment formats (PAF, PSL, trace points).

use crate::error::{FastGAError, Result};
use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

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

    /// Optional PAF tags (everything after column 12)
    pub tags: Vec<String>,

    /// Number of mismatches
    pub mismatches: usize,

    /// Number of gap opens
    pub gap_opens: usize,

    /// Total gap length
    pub gap_len: usize,
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
        let mut paf = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
            self.mapping_quality
        );

        // Add all tags (including CIGAR if present)
        if !self.tags.is_empty() {
            paf.push('\t');
            paf.push_str(&self.tags.join("\t"));
        } else if !self.cigar.is_empty() {
            // If no tags but we have CIGAR, add it as a tag
            paf.push_str(&format!("\tcg:Z:{}", self.cigar));
        }

        paf
    }

    /// Parses a PAF line into an Alignment.
    pub fn from_paf_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            return Err(FastGAError::PafParseError(format!(
                "PAF line has {} fields, expected at least 12",
                fields.len()
            )));
        }

        let mut alignment = Alignment {
            query_name: fields[0].to_string(),
            query_len: fields[1]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid query length".to_string()))?,
            query_start: fields[2]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid query start".to_string()))?,
            query_end: fields[3]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid query end".to_string()))?,
            strand: fields[4]
                .chars()
                .next()
                .ok_or_else(|| FastGAError::PafParseError("Invalid strand".to_string()))?,
            target_name: fields[5].to_string(),
            target_len: fields[6]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid target length".to_string()))?,
            target_start: fields[7]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid target start".to_string()))?,
            target_end: fields[8]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid target end".to_string()))?,
            matches: fields[9]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid matches".to_string()))?,
            block_len: fields[10]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid block length".to_string()))?,
            mapping_quality: fields[11]
                .parse()
                .map_err(|_| FastGAError::PafParseError("Invalid mapping quality".to_string()))?,
            cigar: String::new(),
            mismatches: 0,
            gap_opens: 0,
            gap_len: 0,
            tags: Vec::new(),
        };

        // Parse optional tags and preserve them all
        for field in &fields[12..] {
            if let Some(cigar) = field.strip_prefix("cg:Z:") {
                // Extract CIGAR string
                alignment.cigar = cigar.to_string();
                alignment.parse_cigar_stats()?;
            }
            // Always preserve the tag as-is
            alignment.tags.push(field.to_string());
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
                let count: usize = num_str.parse().map_err(|_| {
                    FastGAError::CigarParseError("Invalid number in CIGAR".to_string())
                })?;
                num_str.clear();

                match ch {
                    'X' => mismatches += count, // Mismatch
                    'I' | 'D' => {
                        // Insertion or Deletion
                        gap_len += count;
                        if !in_gap {
                            gap_opens += 1;
                            in_gap = true;
                        }
                    }
                    '=' | 'M' => in_gap = false, // Match
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
#[derive(Debug, Clone)]
pub struct Alignments {
    /// Vector of individual alignments
    pub alignments: Vec<Alignment>,
}

impl Default for Alignments {
    fn default() -> Self {
        Self::new()
    }
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
        let lines: Vec<String> = self.alignments.iter().map(|a| a.to_paf()).collect();
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
        self.alignments
            .sort_by_key(|a| (a.query_name.clone(), a.query_start));
    }

    /// Sorts alignments by target position.
    pub fn sort_by_target(&mut self) {
        self.alignments
            .sort_by_key(|a| (a.target_name.clone(), a.target_start));
    }

    /// Groups alignments by query sequence.
    pub fn group_by_query(&self) -> std::collections::HashMap<String, Vec<&Alignment>> {
        let mut groups = std::collections::HashMap::new();

        for alignment in &self.alignments {
            groups
                .entry(alignment.query_name.clone())
                .or_insert_with(Vec::new)
                .push(alignment);
        }

        groups
    }

    /// Writes alignments to a file in PAF format.
    ///
    /// # Example
    /// ```no_run
    /// # use fastga_rs::Alignments;
    /// # use anyhow::Result;
    /// # fn main() -> Result<()> {
    /// # let alignments = Alignments::new();
    /// alignments.write_paf("output.paf")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_paf<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        for alignment in &self.alignments {
            writeln!(writer, "{}", alignment.to_paf())?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Writes alignments to a file in tabular format (TSV).
    ///
    /// Includes headers for easy parsing in spreadsheet programs.
    pub fn write_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "query\tq_len\tq_start\tq_end\tstrand\ttarget\tt_len\tt_start\tt_end\tmatches\tblock_len\tidentity\tcigar")?;

        // Write alignments
        for a in &self.alignments {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
                a.query_name,
                a.query_len,
                a.query_start,
                a.query_end,
                a.strand,
                a.target_name,
                a.target_len,
                a.target_start,
                a.target_end,
                a.matches,
                a.block_len,
                a.identity(),
                a.cigar
            )?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Writes alignments to a file in JSON format.
    ///
    /// Useful for programmatic parsing and integration with other tools.
    pub fn write_json<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "[")?;
        for (i, alignment) in self.alignments.iter().enumerate() {
            write!(writer, "  {{")?;
            write!(writer, r#""query_name":"{}","#, alignment.query_name)?;
            write!(writer, r#""query_len":{},"#, alignment.query_len)?;
            write!(writer, r#""query_start":{},"#, alignment.query_start)?;
            write!(writer, r#""query_end":{},"#, alignment.query_end)?;
            write!(writer, r#""strand":"{}","#, alignment.strand)?;
            write!(writer, r#""target_name":"{}","#, alignment.target_name)?;
            write!(writer, r#""target_len":{},"#, alignment.target_len)?;
            write!(writer, r#""target_start":{},"#, alignment.target_start)?;
            write!(writer, r#""target_end":{},"#, alignment.target_end)?;
            write!(writer, r#""matches":{},"#, alignment.matches)?;
            write!(writer, r#""block_len":{},"#, alignment.block_len)?;
            write!(writer, r#""identity":{:.4},"#, alignment.identity())?;
            write!(
                writer,
                r#""mapping_quality":{},"#,
                alignment.mapping_quality
            )?;
            write!(writer, r#""mismatches":{},"#, alignment.mismatches)?;
            write!(writer, r#""gap_opens":{},"#, alignment.gap_opens)?;
            write!(writer, r#""gap_len":{},"#, alignment.gap_len)?;
            write!(writer, r#""cigar":"{}""#, alignment.cigar)?;
            write!(writer, "}}")?;
            if i < self.alignments.len() - 1 {
                writeln!(writer, ",")?;
            } else {
                writeln!(writer)?;
            }
        }
        writeln!(writer, "]")?;

        writer.flush()?;
        Ok(())
    }

    /// Reads alignments from a PAF file.
    pub fn read_paf<P: AsRef<Path>>(path: P) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        Self::from_paf(&contents)
    }

    /// Returns summary statistics about the alignments.
    pub fn summary(&self) -> AlignmentSummary {
        if self.alignments.is_empty() {
            return AlignmentSummary::default();
        }

        let total_matches: usize = self.alignments.iter().map(|a| a.matches).sum();
        let total_mismatches: usize = self.alignments.iter().map(|a| a.mismatches).sum();
        let total_gaps: usize = self.alignments.iter().map(|a| a.gap_len).sum();

        let identities: Vec<f64> = self.alignments.iter().map(|a| a.identity()).collect();
        let mean_identity = identities.iter().sum::<f64>() / identities.len() as f64;
        let min_identity = identities.iter().fold(1.0f64, |a, &b| a.min(b));
        let max_identity = identities.iter().fold(0.0f64, |a, &b| a.max(b));

        AlignmentSummary {
            num_alignments: self.alignments.len(),
            total_matches,
            total_mismatches,
            total_gaps,
            mean_identity,
            min_identity,
            max_identity,
        }
    }
}

/// Summary statistics for a collection of alignments
#[derive(Debug, Clone, Default)]
pub struct AlignmentSummary {
    pub num_alignments: usize,
    pub total_matches: usize,
    pub total_mismatches: usize,
    pub total_gaps: usize,
    pub mean_identity: f64,
    pub min_identity: f64,
    pub max_identity: f64,
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
