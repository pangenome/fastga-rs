#!/bin/bash

# Test FastGA directly
echo "Testing FastGA directly..."

# Create test file
echo -e ">seq1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" > /tmp/test.fa

# Set PATH to include our utilities
export PATH=/home/erik/fastga-rs/target/debug/build/fastga-rs-4b5aee3020c8890b/out:$PATH

# Run FastGA
echo "Running FastGA with PATH set..."
cd /tmp
/home/erik/fastga-rs/target/debug/build/fastga-rs-4b5aee3020c8890b/out/FastGA -pafx test.fa test.fa

echo "Exit code: $?"