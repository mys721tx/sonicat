//! Statistical validation of the sonication and mutation models.
//!
//! These tests check that simulated output follows the documented
//! distributions. Unlike the unit tests, they are probabilistic: each check
//! has a small, known chance of failing even when the code is correct. They
//! are therefore `#[ignore]`d by default and run on their own:
//!
//! ```text
//! cargo test --release --test statistical -- --ignored
//! ```
//!
//! # Error budget
//!
//! * Moment checks (means, variances, rates) accept observations within
//!   [`K_SIGMA`] = 5 standard errors of the expected value, a two-sided
//!   false-alarm rate of ≈ 5.7 × 10⁻⁷ per check under the normal
//!   approximation.
//! * Goodness-of-fit checks use χ² critical values at p = 10⁻⁶.
//!
//! With about 35 checks in total, a correct implementation fails a run with
//! probability ≲ 3 × 10⁻⁵. A failure is therefore strong evidence of a
//! modelling bug; rerun with a different seed to confirm.
//!
//! # Seeds
//!
//! Every test derives its RNG from `SONICAT_STAT_SEED` (default
//! [`DEFAULT_SEED`]) plus a per-test offset, and prints the seed it used.
//! Results are reproducible for a given seed on a given platform and
//! dependency set; the stream may differ across platforms because sampling
//! uses the system math library.

use bio::io::fasta;
use rand::{RngExt, SeedableRng, rngs::StdRng};
use sonicat::{Fragment, Mutator, Sonicator};
use std::io::Write;
use std::process::{Command, Stdio};

/// Number of standard errors a moment may deviate from its expectation.
const K_SIGMA: f64 = 5.0;

/// χ² critical values at p = 10⁻⁶, indexed by degrees of freedom.
const CHI2_CRIT_1E6: [f64; 10] = [
    f64::NAN,
    23.93,
    27.63,
    30.66,
    33.38,
    35.89,
    38.26,
    40.52,
    42.70,
    44.81,
];

const DEFAULT_SEED: u64 = 20_261_117;

/// The RNG for a test, seeded from `SONICAT_STAT_SEED` plus `offset`.
fn rng(offset: u64) -> StdRng {
    let base = std::env::var("SONICAT_STAT_SEED")
        .ok()
        .map(|s| s.parse().expect("SONICAT_STAT_SEED must be a u64"))
        .unwrap_or(DEFAULT_SEED);
    let seed = base.wrapping_add(offset);
    eprintln!("seed = {seed} (SONICAT_STAT_SEED = {base})");
    StdRng::seed_from_u64(seed)
}

fn random_dna(n: usize, rng: &mut StdRng) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
}

#[track_caller]
fn assert_within_sigma(what: &str, observed: f64, expected: f64, se: f64) {
    assert!(
        (observed - expected).abs() <= K_SIGMA * se,
        "{what}: observed {observed}, expected {expected} ± {K_SIGMA}×{se}"
    );
}

/// Pearson's χ² statistic of `observed` counts against `proportions`.
fn chi_squared(observed: &[usize], proportions: &[f64]) -> f64 {
    let n: usize = observed.iter().sum();
    observed
        .iter()
        .zip(proportions)
        .map(|(&o, &p)| {
            let e = n as f64 * p;
            (o as f64 - e).powi(2) / e
        })
        .sum()
}

#[track_caller]
fn assert_fits(what: &str, observed: &[usize], proportions: &[f64]) {
    let df = observed.len() - 1;
    let chi2 = chi_squared(observed, proportions);
    assert!(
        chi2 < CHI2_CRIT_1E6[df],
        "{what}: χ²({df}) = {chi2} ≥ {}, counts = {observed:?}",
        CHI2_CRIT_1E6[df]
    );
}

fn mean(xs: impl IntoIterator<Item = f64>) -> f64 {
    let (sum, n) = xs
        .into_iter()
        .fold((0.0, 0usize), |(s, n), x| (s + x, n + 1));
    sum / n as f64
}

/// Per-base coverage of a sequence of `seq_len` bases by `fragments`.
fn coverage(seq_len: usize, fragments: &[Fragment]) -> Vec<u32> {
    let mut diff = vec![0i64; seq_len + 1];
    for f in fragments {
        diff[f.start] += 1;
        diff[f.end] -= 1;
    }
    let mut acc = 0;
    diff[..seq_len]
        .iter()
        .map(|d| {
            acc += d;
            acc as u32
        })
        .collect()
}

// ============================================================================
// Sonication
// ============================================================================

#[test]
#[ignore = "statistical; run with --ignored"]
fn fragment_count_is_poisson() {
    // N ~ Poisson(λ), λ = depth·L/mean. Over K draws the sample mean has
    // SE √(λ/K) and the sample variance has SE ≈ √((λ + 2λ²)/K).
    let s = Sonicator::new(5.0, 100.0, 10.0).unwrap();
    let seq_len = 2_000;
    let lambda = s.expected_fragments(seq_len);
    let k = 2_000;
    let mut r = rng(1);
    let counts: Vec<f64> = (0..k)
        .map(|_| s.fragments(seq_len, &mut r).len() as f64)
        .collect();
    let m = mean(counts.iter().copied());
    let var = mean(counts.iter().map(|c| (c - m).powi(2)));
    let kf = k as f64;

    assert_within_sigma("mean count", m, lambda, (lambda / kf).sqrt());
    assert_within_sigma(
        "count variance",
        var,
        lambda,
        ((lambda + 2.0 * lambda * lambda) / kf).sqrt(),
    );
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn fragment_length_is_normal() {
    // Lengths ~ round(Normal(μ, σ)); clamping is negligible at μ/σ = 7.5.
    let (mu, sigma) = (300.0, 40.0);
    let s = Sonicator::new(50.0, mu, sigma).unwrap();
    let frags = s.fragments(200_000, &mut rng(2));
    let n = frags.len() as f64;
    let lens: Vec<f64> = frags.iter().map(|f| f.len() as f64).collect();
    let m = mean(lens.iter().copied());
    let sd = mean(lens.iter().map(|l| (l - m).powi(2))).sqrt();

    assert_within_sigma("mean length", m, mu, sigma / n.sqrt());
    // Rounding adds variance 1/12; SE of a normal sample sd is σ/√(2n).
    let expected_sd = (sigma * sigma + 1.0 / 12.0).sqrt();
    assert_within_sigma("length sd", sd, expected_sd, sigma / (2.0 * n).sqrt());
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn total_yield_matches_depth() {
    // Σℓ is compound Poisson: E = λμ = depth·L, Var = λ·E[ℓ²].
    let (depth, mu, sigma, seq_len) = (25.0, 150.0, 15.0, 400_000);
    let s = Sonicator::new(depth, mu, sigma).unwrap();
    let bases: usize = s
        .fragments(seq_len, &mut rng(3))
        .iter()
        .map(Fragment::len)
        .sum();
    let lambda = s.expected_fragments(seq_len);
    let se = (lambda * (mu * mu + sigma * sigma)).sqrt() / seq_len as f64;
    assert_within_sigma("yield / L", bases as f64 / seq_len as f64, depth, se);
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn interior_coverage_matches_depth() {
    for (offset, depth) in [(4, 1.0), (5, 10.0), (6, 50.0)] {
        let (mu, sigma, seq_len, margin) = (150.0, 15.0, 300_000, 500);
        let s = Sonicator::new(depth, mu, sigma).unwrap();
        let cov = coverage(seq_len, &s.fragments(seq_len, &mut rng(offset)));
        let m = mean(cov[margin..seq_len - margin].iter().map(|&c| c as f64));
        // Same compound Poisson argument as the yield, relative to depth.
        let lambda = s.expected_fragments(seq_len);
        let se = depth * ((1.0 + (sigma / mu).powi(2)) / lambda).sqrt();
        assert_within_sigma(&format!("interior coverage at depth {depth}"), m, depth, se);
    }
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn coverage_is_uniform_along_interior() {
    // Each interior block of width B is overlapped by ≈ λ(B + μ)/L fragments.
    let (depth, mu, sigma, seq_len, margin) = (20.0, 100.0, 10.0, 500_000, 500);
    let s = Sonicator::new(depth, mu, sigma).unwrap();
    let cov = coverage(seq_len, &s.fragments(seq_len, &mut rng(7)));
    let interior = &cov[margin..seq_len - margin];
    let width = interior.len() / 10;
    let lambda_block = s.expected_fragments(seq_len) * (width as f64 + mu) / seq_len as f64;
    let se = depth * ((1.0 + (sigma / mu).powi(2)) / lambda_block).sqrt();
    for (i, block) in interior.chunks_exact(width).enumerate() {
        let m = mean(block.iter().map(|&c| c as f64));
        assert_within_sigma(&format!("block {i} coverage"), m, depth, se);
    }
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn start_positions_are_uniform() {
    // With fixed length ℓ, starts are uniform on [0, L-ℓ].
    let (seq_len, len, bins) = (10_150, 150, 10);
    let s = Sonicator::new(200.0, len as f64, 0.0).unwrap();
    let span = seq_len - len + 1;
    let mut counts = vec![0usize; bins];
    for f in s.fragments(seq_len, &mut rng(8)) {
        counts[f.start * bins / span] += 1;
    }
    assert_fits("start bins", &counts, &vec![1.0 / bins as f64; bins]);
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn ends_are_under_covered() {
    // Only fragments starting (ending) exactly at the end cover the first
    // (last) base: Poisson with mean ≈ λ/(L-μ+1) ≈ depth/μ ≈ 0.33 here, so
    // reaching depth/5 = 10 has probability < 10⁻¹².
    let (depth, seq_len) = (50.0, 50_000);
    let s = Sonicator::new(depth, 150.0, 15.0).unwrap();
    let cov = coverage(seq_len, &s.fragments(seq_len, &mut rng(9)));
    for (end, c) in [("first", cov[0]), ("last", cov[seq_len - 1])] {
        assert!((c as f64) < depth / 5.0, "{end} base coverage {c}");
    }
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn sonicat_cli_coverage_matches_depth() {
    let (depth, mu, sigma, len) = (20.0, 150.0, 15.0, 100_000);
    let mut r = rng(10);
    let seq = random_dna(len, &mut r);
    let seed = r.random::<u64>().to_string();

    let mut input = Vec::new();
    fasta::Writer::new(&mut input)
        .write("chr", None, &seq)
        .unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_sonicat"))
        .args(["-d", "20", "-l", "150", "-s", "15", "--seed", &seed])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    // Write stdin concurrently: sonicat emits fragments while still reading,
    // so writing all input first could deadlock once stdout's pipe fills.
    let mut pipe = child.stdin.take().unwrap();
    let writer = std::thread::spawn(move || pipe.write_all(&input));
    let out = child.wait_with_output().unwrap();
    writer.join().unwrap().expect("write sonicat stdin");
    assert!(out.status.success());

    let bases: usize = fasta::Reader::new(&out.stdout[..])
        .records()
        .map(|r| r.unwrap().seq().len())
        .sum();
    let lambda = depth * len as f64 / mu;
    let se = (lambda * (mu * mu + sigma * sigma)).sqrt() / len as f64;
    assert_within_sigma("CLI yield / L", bases as f64 / len as f64, depth, se);
}

// ============================================================================
// Mutation
// ============================================================================

#[derive(Debug, Default)]
struct Events {
    substitution: usize,
    insertion: usize,
    deletion: usize,
    unchanged: usize,
}

/// Classifies the outcome of mutating each base of `seq` individually.
fn classify(m: &Mutator, seq: &[u8], rng: &mut StdRng) -> Events {
    let mut ev = Events::default();
    let mut out = Vec::with_capacity(2);
    for &b in seq {
        out.clear();
        m.mutate_base(b, rng, &mut out);
        match out.as_slice() {
            [] => ev.deletion += 1,
            [x] if *x == b => ev.unchanged += 1,
            [_] => ev.substitution += 1,
            [_, _] => ev.insertion += 1,
            other => panic!("unexpected output {other:?}"),
        }
    }
    ev
}

#[track_caller]
fn assert_rate(what: &str, count: usize, n: usize, p: f64) {
    let nf = n as f64;
    assert_within_sigma(what, count as f64, nf * p, (nf * p * (1.0 - p)).sqrt());
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn event_rates_match_parameters() {
    // Distinct rates so that any mix-up between event types is detected.
    let (s, i, d) = (0.05, 0.10, 0.20);
    let m = Mutator::new(s, i, d).unwrap();
    let n = 200_000;
    let mut r = rng(11);
    let seq = random_dna(n, &mut r);
    let ev = classify(&m, &seq, &mut r);
    assert_rate("substitution", ev.substitution, n, s);
    assert_rate("insertion", ev.insertion, n, i);
    assert_rate("deletion", ev.deletion, n, d);
    assert_rate("unchanged", ev.unchanged, n, 1.0 - s - i - d);
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn default_error_profile_rates() {
    // The Brodin et al. 2013 rates used by `muta`, on a 10 Mb sequence.
    let (s, i, d) = (0.000057, 0.000069, 0.0016);
    let m = Mutator::new(s, i, d).unwrap();
    let n = 10_000_000;
    let mut r = rng(12);
    let seq = random_dna(n, &mut r);
    let ev = classify(&m, &seq, &mut r);
    assert_rate("substitution", ev.substitution, n, s);
    assert_rate("insertion", ev.insertion, n, i);
    assert_rate("deletion", ev.deletion, n, d);
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn substitutions_are_uniform_over_other_nucleotides() {
    let m = Mutator::new(1.0, 0.0, 0.0).unwrap();
    let mut r = rng(13);
    for &from in b"ACGT" {
        let mut counts = [0usize; 4];
        let mut out = Vec::with_capacity(1);
        for _ in 0..30_000 {
            out.clear();
            m.mutate_base(from, &mut r, &mut out);
            counts[b"ACGT".iter().position(|&n| n == out[0]).unwrap()] += 1;
        }
        let others: Vec<usize> = b"ACGT"
            .iter()
            .zip(counts)
            .filter(|&(&n, _)| n != from)
            .map(|(_, c)| c)
            .collect();
        assert_fits(
            &format!("substitutions from {}", from as char),
            &others,
            &[1.0 / 3.0; 3],
        );
    }
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn inserted_nucleotides_are_uniform() {
    let m = Mutator::new(0.0, 1.0, 0.0).unwrap();
    let out = m.mutate(&[b'A'; 100_000], &mut rng(14));
    let mut counts = [0usize; 4];
    for b in out.iter().skip(1).step_by(2) {
        counts[b"ACGT".iter().position(|n| n == b).unwrap()] += 1;
    }
    assert_fits("inserted nucleotides", &counts, &[0.25; 4]);
}

#[test]
#[ignore = "statistical; run with --ignored"]
fn expected_length_change() {
    // Per-base length change is +1 w.p. i and −1 w.p. d, so
    // E[|out|] = n(1 + i − d) and Var = n(i + d − (i − d)²).
    let (s, i, d) = (0.01, 0.03, 0.07);
    let m = Mutator::new(s, i, d).unwrap();
    let n = 500_000;
    let mut r = rng(15);
    let seq = random_dna(n, &mut r);
    let out = m.mutate(&seq, &mut r);
    let nf = n as f64;
    let var = i + d - (i - d).powi(2);
    assert_within_sigma(
        "output length",
        out.len() as f64,
        nf * (1.0 + i - d),
        (nf * var).sqrt(),
    );
}
