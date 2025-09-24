# FastGA-rs Integration Guide for SweepGA

## Overview

This guide explains how to integrate FastGA-rs directly into SweepGA to handle repetitive genomes efficiently without generating massive intermediate PAF files.

## The Problem with Repetitive Genomes

When aligning highly repetitive genomes, traditional workflows suffer from:
- **Disk explosion**: One query can align to thousands of locations, generating 100GB+ PAF files
- **Memory bottleneck**: Loading and processing huge alignment files
- **Unfair distribution**: Some queries may dominate the output while others are underrepresented
- **I/O overhead**: Writing then reading massive intermediate files

## The Solution: Integrated Query-wise Processing

FastGA-rs provides an integrated API that allows SweepGA to:
1. Process queries one-by-one against all targets
2. Apply plane sweep filtering immediately per query
3. Ensure fair mapping distribution
4. Minimize memory and disk usage

## Integration Steps

### 1. Add FastGA-rs as a Dependency

In SweepGA's `Cargo.toml`:

```toml
[dependencies]
# Option 1: Local development
fastga-rs = { path = "../fastga-rs" }

# Option 2: From GitHub
fastga-rs = { git = "https://github.com/pangenome/fastga-rs.git" }

# Option 3: From crates.io (once published)
fastga-rs = "0.1"
```

### 2. Modify SweepGA's Main Processing Pipeline

In `sweepga/src/paf_filter.rs`:

```rust
use fastga_rs::integrated::IntegratedAligner;
use fastga_rs::plane_sweep::PlaneSweepConfig;
use fastga_rs::Config;

pub struct PafFilter {
    // ... existing fields ...
    aligner: Option<IntegratedAligner>,  // NEW: embedded aligner
}

impl PafFilter {
    /// Create a new filter with integrated FastGA support
    pub fn new_with_aligner(params: FilterParams) -> Self {
        let align_config = Config::builder()
            .min_identity(params.min_identity)
            .min_alignment_length(params.min_length)
            .num_threads(params.threads)
            .build();

        let sweep_config = PlaneSweepConfig {
            max_per_query: params.num_mappings.unwrap_or(100),
            max_per_target: params.max_per_target.unwrap_or(1),
            min_identity: params.min_identity,
            min_length: params.min_length,
            max_overlap: 0.5,
        };

        let aligner = IntegratedAligner::new(align_config)
            .with_plane_sweep(sweep_config);

        PafFilter {
            // ... existing fields ...
            aligner: Some(aligner),
        }
    }

    /// Process using embedded FastGA instead of reading PAF file
    pub fn process_with_fastga(
        &mut self,
        query_file: &Path,
        target_file: &Path,
        output: &mut dyn Write,
    ) -> Result<()> {
        let aligner = self.aligner.as_ref()
            .ok_or("Aligner not configured")?;

        // Process query-by-query with immediate filtering
        aligner.process_all_queries(
            query_file,
            target_file,
            |alignment| {
                // Pre-filter at alignment level
                alignment.identity() >= self.min_identity &&
                alignment.query_end - alignment.query_start >= self.min_length
            },
            |query_name, alignments| {
                // Apply SweepGA's sophisticated filtering per query
                let filtered = self.apply_filters_to_query(query_name, alignments);

                // Write results immediately
                for aln in filtered {
                    writeln!(output, "{}", aln.to_paf()).ok();
                }
            },
        )?;

        Ok(())
    }
}
```

### 3. Update Command-line Interface

In `sweepga/src/main.rs`:

```rust
#[derive(Parser)]
struct Args {
    // ... existing args ...

    /// Use integrated FastGA aligner instead of reading PAF
    #[arg(long)]
    use_fastga: bool,

    /// Query FASTA file (when using integrated aligner)
    #[arg(long, required_if("use_fastga", "true"))]
    query: Option<PathBuf>,

    /// Target FASTA file (when using integrated aligner)
    #[arg(long, required_if("use_fastga", "true"))]
    target: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let mut filter = PafFilter::new_with_aligner(params);

    if args.use_fastga {
        // Direct alignment + filtering, no intermediate PAF
        filter.process_with_fastga(
            &args.query.unwrap(),
            &args.target.unwrap(),
            &mut output,
        )?;
    } else {
        // Traditional PAF input processing
        filter.process_paf(&args.input, &mut output)?;
    }

    Ok(())
}
```

## Usage Examples

### Traditional Workflow (Memory Intensive)
```bash
# Step 1: Generate massive PAF file
wfmash genome1.fa genome2.fa > huge.paf  # 100GB for repetitive genome!

# Step 2: Filter with SweepGA
sweepga -i huge.paf -o filtered.paf       # Finally reduce to 1GB
```

### New Integrated Workflow (Memory Efficient)
```bash
# Single step: Align and filter simultaneously
sweepga --use-fastga --query genome1.fa --target genome2.fa -o filtered.paf
# Only writes the 1GB filtered output!
```

### For Highly Repetitive Genomes
```bash
# Aggressive filtering for extreme cases
sweepga --use-fastga \
    --query repetitive_contigs.fa \
    --target reference.fa \
    --num-mappings 50 \          # Keep only top 50 per query
    --min-identity 0.95 \        # High quality only
    --min-length 10000 \         # Long alignments only
    -o filtered.paf
```

## Advanced Integration: Query-wise Processing

For maximum control over memory usage with huge genomes:

```rust
// In sweepga/src/streaming_processor.rs

use fastga_rs::integrated::IntegratedAligner;

pub struct StreamingProcessor {
    aligner: IntegratedAligner,
    target_aggregator: HashMap<String, Vec<Alignment>>,
}

impl StreamingProcessor {
    /// Process each query contig independently
    pub fn process_incrementally(
        &mut self,
        query_file: &Path,
        target_file: &Path,
    ) -> Result<()> {
        // Read queries one at a time
        for (query_name, query_seq) in read_fasta_lazy(query_file) {
            println!("Processing query: {} ({} bp)", query_name, query_seq.len());

            // Align this single query against all targets
            let alignments = self.aligner.align_query_contig(
                &query_name,
                &query_seq,
                target_file,
                |aln| aln.identity() > 0.9,  // Pre-filter
            )?;

            // Apply plane sweep immediately for this query
            let filtered = self.apply_query_plane_sweep(alignments);

            // Aggregate by target for later processing
            for aln in filtered {
                self.target_aggregator
                    .entry(aln.target_name.clone())
                    .or_default()
                    .push(aln);
            }

            // Optionally flush to disk if memory grows too large
            if self.should_flush() {
                self.flush_aggregated_alignments()?;
            }
        }

        // Final target-side plane sweep
        self.apply_target_plane_sweep()?;

        Ok(())
    }
}
```

## Memory and Performance Considerations

### Memory Usage Comparison

| Scenario | Traditional (PAF-based) | Integrated (FastGA-rs) |
|----------|------------------------|------------------------|
| 1GB repetitive genome | ~100GB intermediate PAF | <1GB peak memory |
| Processing model | Load all → Filter | Stream → Filter → Write |
| Disk I/O | Write 100GB → Read 100GB | Write 1GB only |
| Query fairness | First-come-first-served | Guaranteed equal opportunity |

### Performance Tips

1. **Adjust `max_per_query`** based on genome repetitiveness:
   - Low repetition: 100-500 alignments per query
   - High repetition: 10-50 alignments per query
   - Extreme repetition: 5-10 alignments per query

2. **Use `min_identity` to pre-filter**:
   - Set to 0.90+ for high-quality assemblies
   - Set to 0.85+ for draft assemblies
   - Reduces processing overhead early

3. **Enable multi-threading**:
   - FastGA supports parallel alignment
   - Set `num_threads` to available cores

## API Reference

### IntegratedAligner

The main integration point for SweepGA:

```rust
pub struct IntegratedAligner {
    align_config: Config,
    sweep_config: Option<PlaneSweepConfig>,
}

impl IntegratedAligner {
    /// Create a new integrated aligner
    pub fn new(align_config: Config) -> Self;

    /// Enable plane sweep pre-filtering
    pub fn with_plane_sweep(self, config: PlaneSweepConfig) -> Self;

    /// Process a single query contig
    pub fn align_query_contig<F>(
        &self,
        query_name: &str,
        query_seq: &[u8],
        target_file: &Path,
        callback: F,
    ) -> Result<Vec<Alignment>>
    where F: FnMut(&Alignment) -> bool;

    /// Process all queries from a FASTA file
    pub fn process_all_queries<F, G>(
        &self,
        query_file: &Path,
        target_file: &Path,
        per_alignment_callback: F,
        per_query_callback: G,
    ) -> Result<usize>
    where
        F: FnMut(&Alignment) -> bool,
        G: FnMut(&str, Vec<Alignment>);
}
```

### PlaneSweepConfig

Configuration for immediate filtering:

```rust
pub struct PlaneSweepConfig {
    /// Maximum alignments to keep per query
    pub max_per_query: usize,

    /// Maximum alignments per query-target pair
    pub max_per_target: usize,

    /// Minimum identity threshold
    pub min_identity: f64,

    /// Minimum alignment length
    pub min_length: usize,

    /// Maximum overlap fraction (0.0-1.0)
    pub max_overlap: f64,
}
```

## Benefits Summary

1. **No Intermediate Files**: Direct alignment → filter → output pipeline
2. **Memory Efficient**: Process queries incrementally, never load all alignments
3. **Fair Distribution**: Each query gets equal opportunity to map
4. **Embedded Binary**: No external dependencies, single Rust executable
5. **Flexible Filtering**: Apply both plane sweep and SweepGA's advanced filters
6. **Streaming Support**: Can process genomes larger than available RAM

## Next Steps

1. Clone FastGA-rs: `git clone https://github.com/pangenome/fastga-rs.git`
2. Add as dependency to SweepGA
3. Implement the integration points shown above
4. Test with repetitive genomes to verify memory savings
5. Benchmark against traditional PAF-based workflow

## Support

For issues or questions about integration:
- FastGA-rs: https://github.com/pangenome/fastga-rs/issues
- SweepGA: https://github.com/ekg/sweepga/issues