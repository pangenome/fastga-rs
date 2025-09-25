/// Build script for compiling FastGA as a static library.
///
/// This compiles FastGA's C code and links it directly into our Rust binary.
use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=deps/fastga");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let fastga_dir = manifest_dir.join("deps").join("fastga");

    // Compile FastGA as a static library
    println!("cargo:warning=Building FastGA as static library...");

    // We need to compile each program separately to avoid redefining main multiple times
    // First compile FastGA with main renamed
    let mut fastga_build = cc::Build::new();
    fastga_build
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
    let mut fatogdb_build = cc::Build::new();
    fatogdb_build
        .file(fastga_dir.join("FAtoGDB.c"))
        .define("main", "fatogdb_main")
        .include(&fastga_dir)
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("_GNU_SOURCE", None)
        .compile("fatogdb_main");

    // Compile GIXmake with main renamed
    let mut gixmake_build = cc::Build::new();
    gixmake_build
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
    let mut common_build = cc::Build::new();
    common_build
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

    // Also build the utility programs as separate binaries
    // FastGA's system() calls need these
    use std::process::Command;

    let utilities = ["FastGA", "FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF"];

    println!("cargo:warning=Building FastGA utilities...");

    for utility in &utilities {
        let status = Command::new("make")
            .current_dir(&fastga_dir)
            .arg(utility)
            .status()
            .unwrap_or_else(|_| panic!("Failed to build {utility}"));

        if !status.success() {
            panic!("Failed to build {utility}");
        }

        // Copy to OUT_DIR so FastGA can find them
        let src = fastga_dir.join(utility);
        let dst = out_dir.join(utility);
        if src.exists() {
            std::fs::copy(&src, &dst).expect(&format!("Failed to copy {} to OUT_DIR", utility));

            // Make executable
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                let mut perms = std::fs::metadata(&dst).unwrap().permissions();
                perms.set_mode(0o755);
                std::fs::set_permissions(&dst, perms).unwrap();
            }
        }
    }

    println!("cargo:rustc-env=FASTGA_BIN_DIR={}", out_dir.display());
}
