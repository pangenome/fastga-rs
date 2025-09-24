/// Build script for embedding FastGA directly into our Rust binary.
///
/// This compiles FastGA and embeds the binaries as static resources,
/// so users get a single self-contained Rust binary with no external dependencies.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=deps/fastga");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    // Step 1: Build FastGA binaries
    println!("Building FastGA binaries...");

    let fastga_dir = manifest_dir.join("deps").join("fastga");

    // Run make in FastGA directory
    let make_status = Command::new("make")
        .current_dir(&fastga_dir)
        .args(&["FastGA", "ALNtoPAF", "FAtoGDB", "-j4"])
        .status()
        .expect("Failed to build FastGA");

    if !make_status.success() {
        panic!("Failed to compile FastGA binaries");
    }

    // Step 2: Copy binaries to OUT_DIR
    let binaries = ["FastGA", "ALNtoPAF", "FAtoGDB"];
    for binary in &binaries {
        let src = fastga_dir.join(binary);
        let dst = out_dir.join(binary);

        if src.exists() {
            std::fs::copy(&src, &dst)
                .expect(&format!("Failed to copy {} binary", binary));
            println!("cargo:warning=Copied {} to OUT_DIR", binary);
        }
    }

    // Step 3: Tell Rust to embed these binaries
    println!("cargo:rustc-env=FASTGA_BINARY={}", out_dir.join("FastGA").display());
    println!("cargo:rustc-env=ALNTOPAF_BINARY={}", out_dir.join("ALNtoPAF").display());
    println!("cargo:rustc-env=FATOGDB_BINARY={}", out_dir.join("FAtoGDB").display());

    // Step 4: Also compile FastGA as a library (optional, for FFI approach)
    let core_sources = vec![
        "deps/fastga/FastGA.c",
        "deps/fastga/align.c",
        "deps/fastga/alncode.c",
        "deps/fastga/GDB.c",
        "deps/fastga/gene_core.c",
        "deps/fastga/libfastk.c",
        "deps/fastga/ONElib.c",
        "deps/fastga/RSDsort.c",
        "deps/fastga/MSDsort.c",
        "deps/fastga/hash.c",
        "deps/fastga/select.c",
        "deps/fastga/ALNtoPAF.c",
    ];

    let mut build = cc::Build::new();
    build
        .files(&core_sources)
        .include("deps/fastga")
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .define("HAVE_PTHREAD", None);

    // Required system libraries
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=z");

    build.compile("fastga");

    println!("cargo:rustc-link-search=native={}", out_dir.display());
}