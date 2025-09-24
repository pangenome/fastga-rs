# Maintenance Guide

## Updating FastGA

FastGA is included as a git subtree (not submodule) for crate compatibility.

To update FastGA to the latest version:

```bash
# Pull latest changes from FastGA (squashed to keep history clean)
git subtree pull --prefix=deps/fastga https://github.com/thegenemyers/FASTGA.git main --squash

# Remove large example files if they were re-added
rm -rf deps/fastga/EXAMPLE deps/fastga/*.png deps/fastga/*.gz

# Commit the cleanup
git add -A && git commit -m "Clean up after FastGA update"
```

## Why Subtree Instead of Submodule?

- **Crate compatibility**: When published to crates.io or used as a git dependency, subtrees work seamlessly while submodules don't
- **Single repository**: Users don't need to worry about submodule initialization
- **Easier CI/CD**: No extra steps needed in build pipelines
- **Squashed history**: Keeps our git history clean while still allowing updates

## Repository Maintenance

To keep the repository size manageable:

1. Always remove example files and test data from deps/fastga
2. Use `--squash` when pulling subtree updates
3. Periodically run: `git gc --prune=now --aggressive`

## Building FastGA

The build.rs script automatically compiles FastGA when building the crate:
- Detects changed source files
- Compiles only what's needed
- Embeds binaries into the Rust library