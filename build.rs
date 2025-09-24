/// Build script for compiling FastGA C code and linking it with our Rust wrapper.
///
/// This script compiles the necessary FastGA C source files into a static library
/// that can be linked with our Rust code. We only compile the core components
/// needed for alignment functionality, not all the utility programs.

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=deps/fastga");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Core FastGA source files needed for alignment functionality
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
        "deps/fastga_lib.c",  // Our streaming interface
    ];

    // Build the FastGA library
    let mut build = cc::Build::new();

    build
        .files(&core_sources)
        .include("deps/fastga")
        .flag("-O3")
        .flag("-fno-strict-aliasing")
        .warnings(false) // FastGA has some warnings we'll ignore
        .define("HAVE_PTHREAD", None);

    // Add required system libraries flags
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=z");

    // Compile the library
    build.compile("fastga");

    println!("cargo:rustc-link-search=native={}", out_dir.display());
}