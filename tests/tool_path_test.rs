//! GIXmake runs FAtoGDB through system(); it must find it without PATH help.
use fastga_rs::orchestrator::FastGAOrchestrator;

#[test]
fn create_index_finds_fatogdb_without_path() {
    let dir = tempfile::tempdir().unwrap();
    let mut seq = String::new();
    let mut x: u64 = 1;
    for _ in 0..5000 {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seq.push(b"ACGT"[(x >> 62) as usize] as char);
    }
    let fasta = dir.path().join("x.fa");
    std::fs::write(&fasta, format!(">a\n{seq}\n>b\n{seq}\n")).unwrap();
    // Only this test runs in this binary, so changing PATH is safe.
    std::env::set_var("PATH", "/usr/bin:/bin");
    let orchestrator = FastGAOrchestrator {
        temp_dir: dir.path().to_string_lossy().into_owned(),
        ..FastGAOrchestrator::new()
    };
    orchestrator
        .create_index(fasta.to_str().unwrap(), 10)
        .expect("GIXmake must find FAtoGDB");
    assert!(dir.path().join("x.gix").exists(), "no x.gix written");
}
