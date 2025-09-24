# FASTGA Quick Start Guide

## Base Case Alignment: aEleCoq1 vs rDibSmi1

This documents the standard workflow for running FASTGA alignments between two genomes.

### Prerequisites
1. Build FASTGA: `make clean && make -j8`
2. Ensure FASTGA tools are in PATH: `export PATH=$PATH:$(pwd)`

### Input Genomes
- **Target (reference):** aEleCoq1.fa.gz (967MB) - Elephantulus edwardii (Cape elephant shrew)
- **Query:** rDibSmi1.fa.gz (712MB) - Dibamus smithi (Smith's blind skink)

### Basic Alignment (PAF without CIGAR)

```bash
# Run FASTGA with default PAF output
PATH=$PATH:$(pwd) ./FastGA -v -T8 aEleCoq1.fa.gz rDibSmi1.fa.gz > aEleCoq1_vs_rDibSmi1.paf 2> fastga.log

# Output: aEleCoq1_vs_rDibSmi1.paf (2.0MB)
# Runtime: ~5 minutes
# Alignments found: 15,458 with average length 729bp
```

### Alignment with CIGAR Strings (PAF with detailed alignment)

```bash
# Run FASTGA with CIGAR strings using X/= format
PATH=$PATH:$(pwd) ./FastGA -v -T8 -pafx aEleCoq1.fa.gz rDibSmi1.fa.gz > aEleCoq1_vs_rDibSmi1_cigar.paf 2> fastga_cigar.log

# Output: aEleCoq1_vs_rDibSmi1_cigar.paf (8.8MB)
# Runtime: ~5 minutes
# Includes: CIGAR strings with = (match), X (mismatch), I (insertion), D (deletion)
```

## Key Options

- `-v`: Verbose mode (progress to stderr)
- `-T8`: Use 8 threads
- `-pafx`: PAF with CIGAR using X/= operators
- `-pafm`: PAF with CIGAR using M operator
- `-pafS`: PAF with long CS strings
- `-pafs`: PAF with short CS strings
- `-psl`: PSL format output
- `-1:output.1aln`: Save in efficient binary format

## Performance Notes

- Index building: ~1 minute per genome (one-time cost with `-k` option to keep)
- Alignment search: ~3-4 minutes for these genome sizes
- CIGAR generation adds minimal overhead
- Memory usage: ~3-4GB peak

## Output Fields

Standard PAF plus:
- `dv:f`: Divergence fraction (e.g., dv:f:.0629 = 6.29% divergent)
- `df:i`: Number of differences in optimal alignment
- `cg:Z`: CIGAR string (when using -pafx or -pafm options)