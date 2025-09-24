# FastGA-WFMash Integration: External Seed Input

## Objective
Add a single command-line flag to WFMash that allows it to consume pre-computed mapping seeds in PAF format, bypassing its internal index building and MinHash seeding phases entirely, while preserving all chaining, filtering, and output functionality.

## Core Concept
WFMash will accept external mapping seeds (from FastGA or any PAF-producing mapper) and process them through its existing chaining/filtering pipeline exactly as if they were generated internally. These seeds become the input to the chaining process, replacing the MinHash-based seed generation.

## Implementation Requirements

### 1. New Command-Line Flag
Add ONE new flag to WFMash:
```
--input-seeds <FILE>
```
- **Purpose**: Read pre-computed mapping seeds from PAF file to seed the mapping process
- **Alias**: Could also be called `--external-seeds`
- **Input**: Standard PAF format file path (or "-" for stdin)
- **Behavior**: When present, WFMash MUST:
  1. Skip ALL index building
  2. Skip ALL MinHash sketching
  3. Skip ALL mapping computation
  4. Read the PAF file
  5. Continue normally with chaining/filtering

### 2. Expected PAF Input Format
The input PAF must contain AT MINIMUM these standard PAF fields:
```
1.  Query sequence name
2.  Query sequence length
3.  Query start (0-based)
4.  Query end
5.  Strand (+/-)
6.  Target sequence name
7.  Target sequence length
8.  Target start (0-based)
9.  Target end
10. Number of matching bases (can be dummy value)
11. Alignment block length
12. Mapping quality (can be 255 for unknown)
```

Optional but useful tags:
- `dv:f:` - Divergence fraction (FastGA provides this)
- `id:f:` - Identity fraction
- `cg:Z:` - CIGAR string

### 3. Required Code Changes

#### 3.1 Command-Line Parsing
In `src/interface/parse_args.hpp`:
```cpp
// Add new parameter to skch::Parameters
struct Parameters {
    // ... existing fields ...
    bool use_external_seeds = false;
    std::string external_seeds_file = "";
};

// Add argument parsing
args::ValueFlag<std::string> input_seeds(parser, "FILE",
    "Use external PAF seeds to seed the mapping process instead of MinHash",
    {"input-seeds"});

if (input_seeds) {
    map_parameters.use_external_seeds = true;
    map_parameters.external_seeds_file = args::get(input_seeds);
    // IMPORTANT: This should skip index building entirely
}
```

#### 3.2 Main Execution Flow
In `src/interface/main.cpp`:
```cpp
if (map_parameters.use_external_seeds) {
    // Path 1: External seeds
    // - Do NOT build index
    // - Do NOT create sketches
    // - Create minimal SequenceIdManager from FASTAs (for name->ID mapping)
    // - Load PAF seeds
    // - Run chaining/filtering on seeds
    // - Output results
} else {
    // Path 2: Normal WFMash operation (existing code)
    // - Build index
    // - Generate MinHash seeds
    // - Run chaining/filtering
    // - Output results
}
```

#### 3.3 PAF Seed Parser Implementation
Create new file `src/map/include/externalSeeder.hpp`:
```cpp
namespace skch {

class ExternalSeeder {
public:
    // Main entry point for processing external seeds
    static void processExternalSeeds(
        const Parameters& param,
        const std::string& seed_file,
        const SequenceIdManager& id_manager,
        std::ostream& output_stream) {

        // 1. Read PAF seed file
        std::vector<MappingResult> all_seeds = loadPAFSeeds(seed_file, id_manager);

        // 2. Group seeds by query sequence
        auto grouped = groupByQuery(all_seeds);

        // 3. For each query, run through EXISTING WFMash chaining/filtering:
        for (auto& [query_id, seeds] : grouped) {
            // Seeds replace MinHash mappings, now process them:

            // a. Apply chaining (if enabled)
            if (param.chain_mappings) {
                chainAndMergeMappings(seeds, param);
            }

            // b. Apply filtering (if enabled)
            if (param.filterMode != filter::NONE) {
                FilterUtils::filterByScaffolds(seeds, ...);
                FilterUtils::filterByGroup(seeds, ...);
                // ... other existing filters ...
            }

            // c. Output using existing output handler
            OutputHandler::reportReadMappings(seeds, ...);
        }
    }

private:
    // Convert single PAF line to MappingResult
    static MappingResult parsePAFLine(
        const std::string& line,
        const SequenceIdManager& id_manager) {

        std::vector<std::string> fields = split(line, '\t');

        MappingResult result;

        // Parse required PAF fields (indices 0-11)
        std::string query_name = fields[0];
        result.queryStartPos = std::stoi(fields[2]);
        result.queryEndPos = std::stoi(fields[3]);
        result.setStrand(fields[4] == "+" ? strnd::FWD : strnd::REV);

        std::string target_name = fields[5];
        result.refStartPos = std::stoi(fields[7]);
        result.refEndPos = std::stoi(fields[8]);

        // Calculate block length
        result.blockLength = result.refEndPos - result.refStartPos;

        // Get sequence IDs from name->ID mapping
        result.refSeqId = id_manager.getTargetId(target_name);
        // Note: querySeqId handled at grouping stage

        // Parse optional tags (field 12 onwards)
        for (size_t i = 12; i < fields.size(); i++) {
            if (fields[i].substr(0, 5) == "dv:f:") {
                // Divergence to identity conversion
                float divergence = std::stof(fields[i].substr(5));
                float identity = 1.0 - divergence;
                result.nucIdentity = static_cast<uint16_t>(identity * 10000);
            }
            else if (fields[i].substr(0, 5) == "id:f:") {
                // Direct identity
                float identity = std::stof(fields[i].substr(5));
                result.nucIdentity = static_cast<uint16_t>(identity * 10000);
            }
        }

        // Set defaults for missing values
        if (result.nucIdentity == 0) {
            result.nucIdentity = 9000; // Default 90% if not specified
        }
        result.kmerComplexity = 100; // Default full complexity
        result.conservedSketches = 0; // Not applicable for external mappings

        return result;
    }
};

}
```

### 4. Testing Instructions

#### Basic Test Command
```bash
# Step 1: Generate mapping seeds with FastGA
FastGA -T8 -pafx target.fa query.fa > seeds.paf

# Step 2: Process seeds through WFMash chaining/filtering
wfmash --input-seeds seeds.paf target.fa query.fa > output.paf

# Compare with normal WFMash output
wfmash target.fa query.fa > normal_output.paf
```

#### What Should Happen
1. WFMash should NOT print any indexing messages
2. WFMash should NOT print any MinHash sketching messages
3. WFMash SHOULD print chaining/filtering messages (if verbose)
4. Output should contain chained/filtered mappings with chain IDs
5. Output format should be identical to normal WFMash output

### 5. Critical Implementation Notes

#### MUST DO:
- Skip ALL index building when `--input-seeds` is present
- Create SequenceIdManager without building index (just parse FASTA headers)
- Preserve ALL existing chaining/filtering logic unchanged
- Support stdin input with `--input-seeds -`
- Handle both single sequence and multi-sequence FASTAs

#### MUST NOT DO:
- Do NOT modify existing chaining algorithms
- Do NOT modify existing filtering logic
- Do NOT change output format
- Do NOT require index files to exist

### 6. Error Handling

The implementation should handle:
- Missing PAF file → Clear error message
- Malformed PAF lines → Skip with warning, continue processing
- Unknown sequence names → Error with sequence name
- Empty input → Graceful exit with message

### 7. Success Criteria

The implementation is successful when:
1. `wfmash --input-seeds <paf> target.fa query.fa` runs without building index
2. Output contains properly chained mappings
3. Chain IDs are present in output (`ch:Z:` tags)
4. All existing WFMash filtering options work with external seeds
5. Performance is faster than full WFMash run (no indexing overhead)

### 8. Future Extensions (NOT required for initial implementation)

- Support for streaming PAF from stdin
- Support for other mapping formats (SAM, PSL)
- Parallel PAF parsing for large files
- Direct pipe from FastGA without temporary file

## Example Integration Workflow

```bash
#!/bin/bash
# Complete FastGA + WFMash pipeline

# Input files
TARGET="genome1.fa"
QUERY="genome2.fa"
THREADS=16

# Step 1: Generate mapping seeds with FastGA
echo "Running FastGA to generate mapping seeds..."
FastGA -T${THREADS} -pafx ${TARGET} ${QUERY} > mapping_seeds.paf

# Step 2: Chain and filter seeds with WFMash
echo "Running WFMash chaining and filtering on seeds..."
wfmash --input-seeds mapping_seeds.paf \
       -t ${THREADS} \
       -c 2000 \     # chain gap
       -n 5 \        # mappings per segment
       -p 85 \       # filter by identity
       ${TARGET} ${QUERY} > final_output.paf

echo "Pipeline complete!"
```

## Key Distinction from Existing Features

**This is NOT the same as `--align-paf`/`--input-paf`:**
- `--align-paf`: Takes already chained mappings and performs alignment (wflign phase)
- `--input-seeds`: Takes raw mapping seeds and performs chaining/filtering (replaces MinHash seeding)

The new flag operates EARLIER in the pipeline:
```
Normal WFMash:    Index → MinHash seeds → Chain/Filter → (optional) Align
With --input-seeds:         External seeds → Chain/Filter → (optional) Align
With --align-paf:                             Chained PAF → Align
```

## Summary

If any part of this specification is unclear:
1. The flag should skip ALL index building and MinHash seeding
2. The flag should preserve ALL chaining/filtering logic
3. Input is standard PAF format representing mapping seeds
4. Output is standard WFMash output with chains

This is a "seed the mapping process from external file" feature, replacing MinHash seed generation with external seeds.