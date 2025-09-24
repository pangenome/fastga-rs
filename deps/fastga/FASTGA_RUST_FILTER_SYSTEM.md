# SweepGA - Fast Genome Alignment with Sophisticated Filtering

## System Overview

**SweepGA** is a new, standalone system that combines FASTGA's alignment capabilities with a Rust reimplementation of wfmash's filtering algorithms. It "sweeps" through raw alignments to produce clean, filtered results. It is completely separate from wfmash and operates as an independent tool.

## What This System Is

**A two-stage pipeline:**
1. **Stage 1**: FASTGA performs fast genome alignment
2. **Stage 2**: Rust library applies sophisticated filtering to FASTGA's output

**Key characteristics:**
- Standalone executable/library (not integrated into wfmash)
- Rust-based filtering implementation (ported from wfmash's C++ algorithms)
- Memory-efficient (32-byte mapping structures)
- Direct replacement for "FASTGA + post-processing scripts" workflows

## How It Works

### Input Flow
```
Genome A (FASTA) ─┐
                  ├─> FASTGA ─> PAF stream ─> SweepGA Filter ─> Filtered PAF
Genome B (FASTA) ─┘
```

### Data Transformation
1. **FASTGA output**: Standard PAF format with alignments
2. **Internal representation**: Compact 32-byte structures (matching wfmash's layout)
3. **Filtered output**: PAF with sophisticated filtering applied

### Filtering Operations Applied
The Rust library reimplements wfmash's filtering pipeline:

- **Chain merging**: Combines nearby alignments into larger blocks
- **Weak mapping removal**: Filters out short/poor quality alignments
- **Group filtering**: Keeps best N alignments per chromosome/contig
- **Length mismatch filtering**: Removes suspicious high-identity short matches
- **Sparsification**: Reduces redundant overlapping alignments
- **Scaffold filtering**: Complex filtering considering genome structure
- **One-to-one filtering**: Identifies reciprocal best matches

## Why This Approach

### Problems It Solves
1. **FASTGA's limited filtering** - Only has basic length/identity filters
2. **wfmash integration complexity** - Hard to integrate FASTGA into wfmash's C++ codebase
3. **Memory efficiency** - Need to handle large genomes without excessive RAM
4. **Performance comparison** - Need clean way to compare approaches

### Advantages
- **Independent development** - Can iterate without touching wfmash
- **Memory safe** - Rust prevents segfaults/leaks
- **Easy testing** - Can directly compare with wfmash output
- **Modern tooling** - Cargo, built-in benchmarks, easy parallelization
- **Clean interfaces** - Clear separation between alignment and filtering

## Comparison with Existing Tools

### vs FASTGA alone
- **FASTGA**: Fast alignment, basic filtering
- **This system**: FASTGA's speed + sophisticated filtering

### vs wfmash
- **wfmash**: Integrated sketching/alignment/filtering in C++
- **This system**: FASTGA alignment + Rust filtering, separate components

### vs FASTGA + scripts
- **Script approach**: Ad-hoc filtering, high memory, slow
- **This system**: Systematic filtering, memory-efficient, fast

## Technical Details

### Memory Layout
Uses exact same 32-byte struct as wfmash:
- 6 × uint32 (positions, lengths, counts)
- 1 × uint16 (scaled identity)
- 2 × uint8 (flags, complexity)

### Performance Targets
- Memory: < 4GB for mammalian genomes
- Speed: Filtering should add < 20% to FASTGA runtime
- Accuracy: Match wfmash's filtering quality

### Implementation Language
- **Rust** for filtering (memory safety, performance)
- **FASTGA** unchanged (C with FastK library)
- **Interface**: PAF text format (universal)

## Use Cases

### Primary Use Case
Researchers who want:
- FASTGA's alignment speed
- wfmash-quality filtering
- Simple, standalone tool

### Example Workflow
```bash
# Old way (FASTGA + manual filtering)
FastGA genome1.fa genome2.fa > raw.paf
awk '$11 > 5000' raw.paf > filtered.paf  # Crude filtering

# New way (SweepGA)
sweepga genome1.fa genome2.fa \
  --min-identity 0.85 \
  --chain-gap 10000 \
  --filter-mode onetoone > filtered.paf  # Sophisticated filtering
```

## Development Status

**Current state**: Design phase

**Implementation plan**:
1. Create Rust project structure
2. Implement PAF parser
3. Port filtering algorithms from wfmash
4. Create test suite comparing with wfmash
5. Optimize performance
6. Release as standalone tool

## Repository Structure

This will be developed in a new repository, separate from both FASTGA and wfmash:

```
sweepga/
├── src/               # Rust source code
├── tests/             # Test data and integration tests
├── benches/           # Performance benchmarks
├── examples/          # Usage examples
└── comparison/        # Scripts to compare with wfmash
```

## Summary

**SweepGA** is a new, hybrid system that:
1. Uses FASTGA for fast alignment (unchanged)
2. "Sweeps" through alignments with sophisticated filtering (Rust implementation)
3. Operates as a standalone tool (not integrated)
4. Provides a clean comparison point with wfmash
5. Offers modern, maintainable filtering solution

The name captures the essence: it sweeps through the raw alignment data, cleaning and filtering to produce high-quality results.