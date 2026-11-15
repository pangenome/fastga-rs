/// Build script for compiling FastGA as a static library.
///
/// This compiles FastGA's C code and links it directly into our Rust binary.
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=deps/fastga");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let fastga_dir = manifest_dir.join("deps").join("fastga");

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
    println!("cargo:bindir={}", out_dir.display());
}
