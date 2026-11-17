//! End-to-end tests of the `sonicat` and `muta` binaries.

use bio::io::fasta;
use rand::{RngExt, SeedableRng, rngs::StdRng};
use std::collections::HashMap;
use std::io::{ErrorKind, Write};
use std::process::{Command, Output, Stdio};

fn random_dna(n: usize, seed: u64) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
}

fn to_fasta(records: &[(&str, Option<&str>, &[u8])]) -> Vec<u8> {
    let mut buf = Vec::new();
    {
        let mut w = fasta::Writer::new(&mut buf);
        for (id, desc, seq) in records {
            w.write(id, *desc, seq).unwrap();
        }
        w.flush().unwrap();
    }
    buf
}

/// Runs `bin` with `args`, feeding it `stdin`, and collects its output.
///
/// Stdin is written from a separate thread while stdout and stderr are
/// drained, so a child that writes output before consuming all its input
/// cannot deadlock on a full pipe. A child may also exit without reading
/// stdin at all (e.g. on invalid arguments), so a broken pipe is expected
/// and ignored.
fn run(bin: &str, args: &[&str], stdin: &[u8]) -> Output {
    let mut child = Command::new(bin)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn binary");
    let mut pipe = child.stdin.take().unwrap();
    let input = stdin.to_vec();
    let writer = std::thread::spawn(move || match pipe.write_all(&input) {
        Err(e) if e.kind() == ErrorKind::BrokenPipe => {}
        result => result.expect("write child stdin"),
    });
    let output = child.wait_with_output().unwrap();
    writer.join().unwrap();
    output
}

fn sonicat(args: &[&str], stdin: &[u8]) -> Output {
    run(env!("CARGO_BIN_EXE_sonicat"), args, stdin)
}

fn muta(args: &[&str], stdin: &[u8]) -> Output {
    run(env!("CARGO_BIN_EXE_muta"), args, stdin)
}

fn parse(stdout: &[u8]) -> Vec<fasta::Record> {
    fasta::Reader::new(stdout)
        .records()
        .map(|r| r.expect("valid FASTA output"))
        .collect()
}

/// Parses a `ID:START-END` description (1-based, inclusive).
fn parse_region(desc: &str) -> (&str, usize, usize) {
    let (id, range) = desc.rsplit_once(':').expect("ID:START-END");
    let (start, end) = range.split_once('-').expect("START-END");
    (id, start.parse().unwrap(), end.parse().unwrap())
}

// ---- sonicat -----------------------------------------------------------------

#[test]
fn sonicat_fragments_match_their_source_region() {
    let chr1 = random_dna(5_000, 1);
    let chr2 = random_dna(3_000, 2);
    let input = to_fasta(&[("chr1", None, &chr1), ("chr2", Some("second"), &chr2)]);
    let out = sonicat(
        &["-d", "10", "-l", "200", "-s", "30", "--seed", "42"],
        &input,
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let sources: HashMap<&str, &[u8]> = [("chr1", &chr1[..]), ("chr2", &chr2[..])].into();
    // Holds for every fragment the run produces, whatever the draw.
    for (n, rec) in parse(&out.stdout).iter().enumerate() {
        assert_eq!(rec.id(), format!("seq_{}", n + 1), "names are sequential");
        let (src, start, end) = parse_region(rec.desc().expect("region description"));
        assert!(start >= 1 && start <= end, "bad region {start}-{end}");
        assert_eq!(rec.seq(), &sources[src][start - 1..end]);
    }
}

#[test]
fn sonicat_seed_is_reproducible() {
    let input = to_fasta(&[("chr", None, &random_dna(2_000, 3))]);
    let args = ["-d", "5", "--seed", "7"];
    let a = sonicat(&args, &input);
    let b = sonicat(&args, &input);
    assert!(a.status.success());
    assert_eq!(a.stdout, b.stdout);
}

#[test]
fn sonicat_streams_input_and_output_larger_than_pipe_buffers() {
    // ~1 MiB in and out: completes only if the harness feeds stdin while
    // draining stdout, since sonicat writes fragments as it reads.
    let input = to_fasta(&[("chr", None, &random_dna(1 << 20, 13))]);
    let out = sonicat(&["-d", "1", "--seed", "1"], &input);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    for rec in parse(&out.stdout) {
        let (_, start, end) = parse_region(rec.desc().unwrap());
        assert!(start >= 1 && start <= end && end <= 1 << 20);
    }
}

#[test]
fn sonicat_reads_and_writes_files() {
    let dir = std::env::temp_dir().join(format!("sonicat-cli-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fin = dir.join("in.fa");
    let fout = dir.join("out.fa");
    let input = to_fasta(&[("chr", None, &random_dna(1_000, 5))]);
    std::fs::write(&fin, &input).unwrap();

    let out = sonicat(
        &[
            "-i",
            fin.to_str().unwrap(),
            "-o",
            fout.to_str().unwrap(),
            "--seed",
            "1",
        ],
        b"",
    );
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
    // Same seed through stdin/stdout must give byte-identical output.
    let piped = sonicat(&["--seed", "1"], &input);
    assert_eq!(std::fs::read(&fout).unwrap(), piped.stdout);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn sonicat_rejects_invalid_parameters() {
    // Larger than any OS pipe buffer, so the child always exits before the
    // harness finishes writing stdin (regression for EPIPE on Linux).
    let input = to_fasta(&[("chr", None, &random_dna(1 << 20, 12))]);
    for args in [
        &["-d", "-1"][..],
        &["-l", "0"],
        &["-s", "-5"],
        &["-d", "abc"],
    ] {
        let out = sonicat(args, &input);
        assert!(!out.status.success(), "accepted {args:?}");
        assert!(!out.stderr.is_empty());
    }
}

#[test]
fn sonicat_empty_input_gives_empty_output() {
    let out = sonicat(&["--seed", "1"], b"");
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
}

// ---- muta --------------------------------------------------------------------

#[test]
fn muta_zero_rates_is_identity() {
    let input = to_fasta(&[
        ("a", Some("first record"), &random_dna(1_000, 6)),
        ("b", None, b"acgtNNNN"),
    ]);
    let out = muta(&["-s", "0", "-n", "0", "-d", "0", "--seed", "1"], &input);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(parse(&out.stdout), parse(&input));
}

#[test]
fn muta_preserves_record_ids_and_descriptions() {
    let input = to_fasta(&[
        ("a", Some("first record"), &random_dna(500, 7)),
        ("b", None, &random_dna(500, 8)),
    ]);
    let out = muta(
        &["-s", "0.1", "-n", "0.1", "-d", "0.1", "--seed", "2"],
        &input,
    );
    assert!(out.status.success());
    let recs = parse(&out.stdout);
    assert_eq!(recs.len(), 2);
    assert_eq!((recs[0].id(), recs[0].desc()), ("a", Some("first record")));
    assert_eq!((recs[1].id(), recs[1].desc()), ("b", None));
}

#[test]
fn muta_certain_deletion_empties_sequences() {
    let input = to_fasta(&[("a", None, &random_dna(500, 9))]);
    let out = muta(&["-s", "0", "-n", "0", "-d", "1", "--seed", "3"], &input);
    assert!(out.status.success());
    let recs = parse(&out.stdout);
    assert_eq!(recs.len(), 1);
    assert!(recs[0].seq().is_empty());
}

#[test]
fn muta_seed_is_reproducible() {
    let input = to_fasta(&[("a", None, &random_dna(5_000, 10))]);
    let args = ["-s", "0.01", "-n", "0.01", "-d", "0.01", "--seed", "4"];
    assert_eq!(muta(&args, &input).stdout, muta(&args, &input).stdout);
}

#[test]
fn muta_rejects_invalid_parameters() {
    // Larger than any OS pipe buffer, so the child always exits before the
    // harness finishes writing stdin (regression for EPIPE on Linux).
    let input = to_fasta(&[("a", None, &random_dna(1 << 20, 12))]);
    for args in [
        &["-s", "-0.1"][..],
        &["-n", "1.5"],
        &["-s", "0.5", "-n", "0.5", "-d", "0.5"],
    ] {
        let out = muta(args, &input);
        assert!(!out.status.success(), "accepted {args:?}");
        assert!(!out.stderr.is_empty());
    }
}
