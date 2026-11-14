use fastga_rs::binary_finder::find_binary;
use fastga_rs::runner::Orchestrator;
use fastga_rs::{AlnReader, AlnWriter};
use std::path::Path;
use std::process::Command;

fn dna(len: usize, mut state: u64) -> Vec<u8> {
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn mutate(seq: &[u8], every: usize) -> Vec<u8> {
    let mut out = seq.to_vec();
    for base in out.iter_mut().step_by(every) {
        *base = if *base == b'A' { b'C' } else { b'A' };
    }
    out
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn write_fasta(path: &Path, records: &[(&str, Vec<u8>)]) {
    let text: String = records
        .iter()
        .map(|(name, seq)| format!(">{name}\n{}\n", String::from_utf8_lossy(seq)))
        .collect();
    std::fs::write(path, text).unwrap();
}

/// PAF lines of a .1aln file, from ALNtoPAF.
fn paf_lines(aln: &Path) -> Vec<String> {
    let bin = find_binary("ALNtoPAF").unwrap();
    let out = Command::new(bin)
        .arg("-T1")
        .arg(aln.with_extension(""))
        .current_dir(aln.parent().unwrap())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "ALNtoPAF: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

#[test]
fn copy_keeps_every_line_of_the_selected_records() {
    let dir = tempfile::tempdir().unwrap();
    let base = dna(6000, 7);
    let query = dir.path().join("query.fa");
    let target = dir.path().join("target.fa");
    write_fasta(&query, &[("q", base.clone())]);
    write_fasta(
        &target,
        &[
            ("t_fwd", mutate(&base, 211)),
            ("t_rev", revcomp(&mutate(&base, 157))),
            ("t_fwd2", mutate(&base, 113)),
        ],
    );
    let produced = Orchestrator::new_simple(1).align(&query, &target).unwrap();
    let input = dir.path().join("input.1aln");
    std::fs::copy(produced.trim(), &input).unwrap();

    // The reader returns every record, the first included
    let read = AlnReader::open(&input).unwrap().read_all().unwrap();
    let total = read.len();
    let input_paf = paf_lines(&input);
    assert!(total >= 3, "need 3 alignments, got {total}");
    assert_eq!(input_paf.len(), total);
    for (aln, line) in read.iter().zip(&input_paf) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(aln.strand.to_string(), f[4]);
        assert_eq!(aln.target_start.to_string(), f[7]);
        assert_eq!(aln.target_end.to_string(), f[8]);
    }
    assert!(input_paf.iter().any(|l| l.split('\t').nth(4) == Some("-")));

    // All records: same PAF as the input
    let all = dir.path().join("all.1aln");
    let mut writer = AlnWriter::create_with_gdb(&all, &input, true).unwrap();
    assert_eq!(
        writer.copy_alignments_from(&input, |_| true).unwrap(),
        total
    );
    writer.finalize();
    assert_eq!(paf_lines(&all), input_paf);

    // One record dropped, last record kept
    let some = dir.path().join("some.1aln");
    let mut writer = AlnWriter::create_with_gdb(&some, &input, true).unwrap();
    assert_eq!(
        writer
            .copy_alignments_from(&input, |rank| rank != 1)
            .unwrap(),
        total - 1
    );
    writer.finalize();
    let mut expected = input_paf.clone();
    expected.remove(1);
    assert_eq!(paf_lines(&some), expected);
}
