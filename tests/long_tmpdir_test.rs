use anyhow::Result;
use fastga_rs::binary_finder::find_binary;
use std::fs::{self, File};
use std::io::Write;
use std::process::Command;
use tempfile::tempdir;

/// ONElib wrote "$TMPDIR/OneSchema.XXXXXX" into a 64-byte buffer, so a $TMPDIR
/// longer than 46 bytes smashed the stack and aborted the process.
#[test]
fn fatogdb_survives_a_long_tmpdir() -> Result<()> {
    let dir = tempdir()?;

    let fasta = dir.path().join("seq.fa");
    let mut file = File::create(&fasta)?;
    writeln!(file, ">sequence1")?;
    writeln!(file, "{}", "ACGTACGTACGTACGTACGTACGTACGTACGT".repeat(40))?;
    file.flush()?;

    // 120 characters of directory name, well past the old 46-byte limit.
    let long_tmp = dir.path().join("t".repeat(120));
    fs::create_dir_all(&long_tmp)?;
    assert!(long_tmp.as_os_str().len() > 46);

    let fatogdb = find_binary("FAtoGDB")?;
    let output = Command::new(&fatogdb)
        .env("TMPDIR", &long_tmp)
        .arg(&fasta)
        .output()?;

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("buffer overflow"),
        "FAtoGDB overflowed its temp-path buffer: {stderr}"
    );
    assert!(
        output.status.success(),
        "FAtoGDB failed with {:?}: {stderr}",
        output.status.code()
    );

    Ok(())
}
