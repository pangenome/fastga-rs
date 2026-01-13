# Updating the FASTGA Git Subtree

## Overview

`fastga-rs` uses a **git subtree** (not submodule) to embed the FASTGA source code at `deps/fastga`. This allows `cargo install` and `cargo publish` to work properly without requiring submodule initialization.

## Current Setup

- **Subtree path**: `deps/fastga`
- **Upstream repository**: `https://github.com/thegenemyers/FASTGA.git`
- **Current commit**: `6c5b8b6` (V1.5, includes ANO annotation support)
- **Branch**: `main`

**Note**: As of 2026-01-13, zstd compression support was removed to match upstream.
Upstream FASTGA removed zstd-compressed ktab support in favor of simpler uncompressed indexes.

## When to Update

Update the FASTGA subtree when:
- Gene Myers releases a new version with bug fixes
- New features are added to FASTGA that you need
- Performance improvements are available

## How to Update

### 1. Check for upstream changes

```bash
# Add the upstream remote (only needed once)
git remote add fastga-upstream https://github.com/thegenemyers/FASTGA.git

# Fetch latest changes
git fetch fastga-upstream main

# View what's new
git log HEAD..fastga-upstream/main --oneline
```

### 2. Update the subtree

```bash
# Pull latest changes from upstream into the subtree
git subtree pull --prefix deps/fastga https://github.com/thegenemyers/FASTGA.git main --squash

# This creates a merge commit with the new FASTGA code
```

### 3. Test the update

```bash
# Clean build to ensure everything compiles
cargo clean
cargo build --release

# Test that GDB skeleton embedding still works
./target/release/build/fastga-rs-*/out/FastGA data/cerevisiae.chrV.fa.gz -1:/tmp/test.1aln
ONEview /tmp/test.1aln | grep -A 5 "^g$"

# Should show GDB skeleton with scaffold names:
# g
# S 12 S288C#1#chrV
# C 583092
# ...
```

### 4. Commit and push

```bash
# The subtree pull creates a commit automatically
# Just push to origin
git push origin main
```

## Updating to a Specific Commit

If you need a specific commit instead of the latest `main`:

```bash
# Update to a specific commit hash
git subtree pull --prefix deps/fastga https://github.com/thegenemyers/FASTGA.git <commit-hash> --squash
```

## Troubleshooting

### "working tree has modifications" error

If you get this error:
```bash
# Commit or stash your changes first
git status
git add -A
git commit -m "Save work before FASTGA update"

# Then retry the subtree pull
git subtree pull --prefix deps/fastga https://github.com/thegenemyers/FASTGA.git main --squash
```

### Build fails after update

If the build fails after updating:
1. Check if FASTGA changed its API or file structure
2. Update `build.rs` if needed to match new file names
3. Check for breaking changes in FASTGA's changelog

### Reverting an update

If the update breaks something:
```bash
# Find the commit before the subtree pull
git log --oneline -5

# Reset to before the update
git reset --hard <commit-before-update>
```

## Alternative: Update from a Fork

If you're using your own FASTGA fork (like `ekg/FASTGA`):

```bash
# Update from your fork instead
git subtree pull --prefix deps/fastga https://github.com/ekg/FASTGA.git <branch> --squash
```

## How This Differs from Submodules

**Git Subtree** (what we use):
- ✅ FASTGA source is committed directly into fastga-rs
- ✅ Works with `cargo install` and `cargo publish`
- ✅ No initialization needed when cloning
- ❌ Slightly more complex to update

**Git Submodule** (old approach):
- ❌ FASTGA is a separate repository reference
- ❌ Doesn't work with `cargo install` or `cargo publish`
- ❌ Requires `git submodule update --init` after cloning
- ✅ Simpler to update

## Verifying GDB Skeleton Support

After any FASTGA update, verify the critical features:

```bash
# Generate a .1aln file
./target/release/build/fastga-rs-*/out/FastGA data/cerevisiae.chrV.fa.gz -1:/tmp/verify.1aln

# Check for GDB skeleton (scaffold names)
ONEview /tmp/verify.1aln | grep -c "^S "
# Should be > 0 (number of scaffolds)

# Check for X fields (edit distance)
ONEview /tmp/verify.1aln | grep "^X " | head -3
# Should show X records with edit distance data
```

## History

- **2026-01-13**: Updated to commit `6c5b8b6` (V1.5), removed zstd support to match upstream
- **2025-10-09**: Converted from git submodule to git subtree
- **2025-10-09**: Updated to commit `a3c1564` for GDB skeleton support
- Previous version used submodule at `ekg/FASTGA` branch `makefile-fixes`

## References

- [Git Subtree Documentation](https://git-scm.com/book/en/v2/Git-Tools-Advanced-Merging)
- [FASTGA Repository](https://github.com/thegenemyers/FASTGA)
- [fastga-rs Issue Tracker](https://github.com/pangenome/fastga-rs/issues)
