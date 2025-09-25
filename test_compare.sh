#!/bin/bash
# Compare our orchestrator vs direct FastGA

# Extract chrV data
gunzip -c data/cerevisiae.chrV.fa.gz > test_chrv.fa

echo "=== Running FastGA directly ==="
# Run FastGA directly (it will call system() for utilities)
./target/debug/build/fastga-rs-4b5aee3020c8890b/out/FastGA \
    -pafx -T4 -l1000 -i0.95 -M \
    test_chrv.fa test_chrv.fa > direct_fastga.paf 2>/dev/null

echo "Direct FastGA: $(wc -l < direct_fastga.paf) alignments"

echo ""
echo "=== Running our Rust orchestrator ==="
# Use our library
cat > test_ours.rs << 'EOF'
use fastga_rs::{FastGA, Config};
use std::path::Path;

fn main() {
    let config = Config::builder()
        .num_threads(4)
        .min_alignment_length(1000)
        .min_identity(0.95)
        .soft_masking(true)
        .build();

    let aligner = FastGA::new(config).unwrap();
    let result = aligner.align_files(
        Path::new("test_chrv.fa"),
        Path::new("test_chrv.fa")
    ).unwrap();

    // Write PAF output
    result.write_paf("our_orchestrator.paf").unwrap();
    println!("Our orchestrator: {} alignments", result.len());
}
EOF

cargo build --example verify_alignment 2>/dev/null
rustc test_ours.rs -L target/debug/deps --extern fastga_rs=target/debug/libfastga_rs.rlib 2>/dev/null && ./test_ours

echo ""
echo "=== Comparing outputs ==="
echo "Direct FastGA alignments: $(wc -l < direct_fastga.paf)"
echo "Our orchestrator alignments: $(wc -l < our_orchestrator.paf)"

# Compare first few lines (should be identical except for runtime differences)
echo ""
echo "First alignment from direct FastGA:"
head -1 direct_fastga.paf | cut -f1-12

echo ""
echo "First alignment from our orchestrator:"
head -1 our_orchestrator.paf | cut -f1-12

# Check if outputs are functionally identical
echo ""
echo "Checking if alignments match..."
# Sort both files and compare (PAF order might vary slightly)
sort direct_fastga.paf > direct_sorted.paf
sort our_orchestrator.paf > our_sorted.paf

if diff -q direct_sorted.paf our_sorted.paf > /dev/null; then
    echo "✅ PERFECT MATCH: Outputs are identical!"
else
    echo "Comparing alignment counts per query:"
    cut -f1 direct_fastga.paf | sort | uniq -c | head -5
    echo "vs"
    cut -f1 our_orchestrator.paf | sort | uniq -c | head -5
fi