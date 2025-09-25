#!/bin/bash

echo "Testing individual utilities..."

BIN_DIR=/home/erik/fastga-rs/target/debug/build/fastga-rs-4b5aee3020c8890b/out
cd /tmp

# Test FAtoGDB
echo -e ">test\nACGTACGTACGTACGT" > test.fa
echo "Running FAtoGDB..."
$BIN_DIR/FAtoGDB test.fa
echo "FAtoGDB exit code: $?"
ls -la test.1gdb 2>/dev/null

# Test GIXmake
echo "Running GIXmake..."
$BIN_DIR/GIXmake -f10 test
echo "GIXmake exit code: $?"
ls -la test.gix 2>/dev/null