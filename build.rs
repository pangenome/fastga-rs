/// Build script for compiling FastGA as a static library.
///
/// This compiles FastGA's C code and links it directly into our Rust binary.
use std::env;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::Command;

/// Patch ONElib.c to respect TMPDIR environment variable instead of hardcoded /tmp/
/// This allows FastGA to work on systems where /tmp is not writable or is on a
/// read-only filesystem (common in HPC environments).
fn patch_onelib_tmpdir(fastga_dir: &std::path::Path) {
    let onelib_path = fastga_dir.join("ONElib.c");

    let mut content = String::new();
    std::fs::File::open(&onelib_path)
        .expect("Failed to open ONElib.c")
        .read_to_string(&mut content)
        .expect("Failed to read ONElib.c");

    // Check if already patched
    if content.contains("get_tmpdir_path") {
        return;
    }

    // Add helper function after the includes section
    // We insert after "static int listEltSize" which appears early in the file
    let helper_function = r#"

// Helper function to get temp directory from TMPDIR env var, fallback to /tmp
// Patched by fastga-rs build script to support HPC environments
static const char* get_tmpdir_path(void) {
    const char* tmpdir = getenv("TMPDIR");
    if (tmpdir && tmpdir[0] != '\0') {
        return tmpdir;
    }
    tmpdir = getenv("TMP");
    if (tmpdir && tmpdir[0] != '\0') {
        return tmpdir;
    }
    tmpdir = getenv("TEMP");
    if (tmpdir && tmpdir[0] != '\0') {
        return tmpdir;
    }
    return "/tmp";
}
"#;

    // Find insertion point after "static int listEltSize" line
    let insertion_marker = "static int listEltSize";
    if let Some(pos) = content.find(insertion_marker) {
        if let Some(newline_pos) = content[pos..].find('\n') {
            let insert_pos = pos + newline_pos + 1;
            content.insert_str(insert_pos, helper_function);
        }
    }

    // Replace hardcoded /tmp/ paths with dynamic temp directory
    // Pattern 1: sprintf (template, "/tmp/OneSchema.%d", getpid()) ;
    content = content.replace(
        r#"sprintf (template, "/tmp/OneSchema.%d", getpid())"#,
        r#"sprintf (template, "%s/OneSchema.%d", get_tmpdir_path(), getpid())"#,
    );

    // Pattern 2: strcpy (template, "/tmp/OneSchema.XXXXXX") ;
    content = content.replace(
        r#"strcpy (template, "/tmp/OneSchema.XXXXXX")"#,
        r#"sprintf (template, "%s/OneSchema.XXXXXX", get_tmpdir_path())"#,
    );

    // Pattern 3: char template[] = "/tmp/OneTextSchema-XXXXXX" ;
    // Need to change fixed-size array to larger buffer for sprintf
    content = content.replace(
        r#"char template[] = "/tmp/OneTextSchema-XXXXXX""#,
        r#"char template[4096]; sprintf(template, "%s/OneTextSchema-XXXXXX", get_tmpdir_path())"#,
    );

    // Write patched content back
    std::fs::File::create(&onelib_path)
        .expect("Failed to create patched ONElib.c")
        .write_all(content.as_bytes())
        .expect("Failed to write patched ONElib.c");

    println!("cargo:warning=Patched ONElib.c to respect TMPDIR environment variable");
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=deps/fastga");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let fastga_dir = manifest_dir.join("deps").join("fastga");

    // Patch ONElib.c to respect TMPDIR before any compilation
    // This must happen before both cc::Build and make commands
    patch_onelib_tmpdir(&fastga_dir);

    println!("cargo:warning=Building FastGA...");

    // Build our wrapper functions
    let deps_dir = manifest_dir.join("deps");
    cc::Build::new()
        .file(deps_dir.join("rust_wrappers.c"))
        .file(deps_dir.join("rust_onelib.c"))
        .file(fastga_dir.join("GDB.c"))
        .file(fastga_dir.join("gene_core.c"))
        .file(fastga_dir.join("ONElib.c"))
        .file(fastga_dir.join("alncode.c"))
        .include(&fastga_dir)
        .include(&manifest_dir)
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("_GNU_SOURCE", None)
        .compile("fastga_wrappers");

    // We need to compile each program separately to avoid redefining main multiple times
    // First compile FastGA with main renamed
    cc::Build::new()
        .file(fastga_dir.join("FastGA.c"))
        .define("main", "fastga_main")
        .file(fastga_dir.join("align.c"))
        .file(fastga_dir.join("alncode.c"))
        .file(fastga_dir.join("RSDsort.c"))
        .include(&fastga_dir)
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("_GNU_SOURCE", None)
        .compile("fastga_main");

    // Compile FAtoGDB with main renamed
    cc::Build::new()
        .file(fastga_dir.join("FAtoGDB.c"))
        .define("main", "fatogdb_main")
        .include(&fastga_dir)
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("_GNU_SOURCE", None)
        .compile("fatogdb_main");

    // Compile GIXmake with main renamed
    cc::Build::new()
        .file(fastga_dir.join("GIXmake.c"))
        .define("main", "gixmake_main")
        .file(fastga_dir.join("MSDsort.c"))
        .include(&fastga_dir)
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("_GNU_SOURCE", None)
        .compile("gixmake_main");

    // Compile common dependencies
    cc::Build::new()
        .file(fastga_dir.join("GDB.c"))
        .file(fastga_dir.join("gene_core.c"))
        .file(fastga_dir.join("libfastk.c"))
        .file(fastga_dir.join("ONElib.c"))
        .file(fastga_dir.join("hash.c"))
        .file(fastga_dir.join("select.c"))
        .include(&fastga_dir)
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("_GNU_SOURCE", None)
        .compile("fastga_common");

    // Link required system libraries
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=z");

    // Build the utility programs as separate binaries
    // FastGA's system() calls need these
    let utilities = [
        "FastGA", "FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF", "PAFtoALN", "ONEview",
    ];

    println!("cargo:warning=Building FastGA utilities...");

    for utility in &utilities {
        let mut cmd = Command::new("make");
        cmd.current_dir(&fastga_dir).arg(utility);

        let status = cmd
            .status()
            .unwrap_or_else(|_| panic!("Failed to build {utility}"));

        if !status.success() {
            panic!("Failed to build {utility}");
        }

        // Copy to OUT_DIR so FastGA can find them
        let src = fastga_dir.join(utility);
        let dst = out_dir.join(utility);
        if src.exists() {
            std::fs::copy(&src, &dst)
                .unwrap_or_else(|_| panic!("Failed to copy {utility} to OUT_DIR"));

            // Make executable
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let mut perms = std::fs::metadata(&dst).unwrap().permissions();
                perms.set_mode(0o755);
                std::fs::set_permissions(&dst, perms).unwrap();
            }

            // Remove from source directory to avoid cargo publish verification errors
            // (cargo doesn't allow build scripts to modify source directory)
            let _ = std::fs::remove_file(&src);
        }
    }

    println!("cargo:rustc-env=FASTGA_BIN_DIR={}", out_dir.display());
}
