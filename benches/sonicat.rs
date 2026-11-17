//! Throughput benchmarks for the sonication and mutation simulators.
//!
//! Run with `cargo bench`; throughput is reported per input base.

use bio::io::fasta;
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::{RngExt, SeedableRng, rngs::StdRng};
use sonicat::{Mutator, Sonicator};
use std::hint::black_box;
use std::io;

fn random_dna(n: usize, seed: u64) -> Vec<u8> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
}

/// Fragment generation alone, across genome sizes at the CLI defaults.
fn sonicate_by_length(c: &mut Criterion) {
    let mut group = c.benchmark_group("sonicate/length");
    let s = Sonicator::new(50.0, 150.0, 15.0).unwrap();
    for len in [10_000, 100_000, 1_000_000] {
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, &len| {
            let mut rng = StdRng::seed_from_u64(1);
            b.iter(|| black_box(s.fragments(black_box(len), &mut rng)));
        });
    }
    group.finish();
}

/// Fragment generation scales with the number of fragments, i.e. depth.
fn sonicate_by_depth(c: &mut Criterion) {
    let mut group = c.benchmark_group("sonicate/depth");
    let len = 100_000;
    group.throughput(Throughput::Bytes(len as u64));
    for depth in [1.0, 10.0, 100.0] {
        let s = Sonicator::new(depth, 150.0, 15.0).unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(depth), &depth, |b, _| {
            let mut rng = StdRng::seed_from_u64(2);
            b.iter(|| black_box(s.fragments(len, &mut rng)));
        });
    }
    group.finish();
}

/// Fixed vs. normally distributed fragment lengths.
fn sonicate_length_model(c: &mut Criterion) {
    let mut group = c.benchmark_group("sonicate/sd");
    let len = 100_000;
    group.throughput(Throughput::Bytes(len as u64));
    for sd in [0.0, 15.0] {
        let s = Sonicator::new(50.0, 150.0, sd).unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(sd), &sd, |b, _| {
            let mut rng = StdRng::seed_from_u64(3);
            b.iter(|| black_box(s.fragments(len, &mut rng)));
        });
    }
    group.finish();
}

/// The full `sonicat` pipeline minus I/O: shear, slice and format as FASTA.
fn sonicate_to_fasta(c: &mut Criterion) {
    let mut group = c.benchmark_group("sonicate/fasta");
    let seq = random_dna(100_000, 4);
    let s = Sonicator::new(50.0, 150.0, 15.0).unwrap();
    group.throughput(Throughput::Bytes(seq.len() as u64));
    group.bench_function("100000", |b| {
        let mut rng = StdRng::seed_from_u64(5);
        b.iter(|| {
            let mut w = fasta::Writer::new(io::sink());
            for (n, f) in s.fragments(seq.len(), &mut rng).iter().enumerate() {
                let name = format!("seq_{}", n + 1);
                let desc = format!("chr:{}-{}", f.start + 1, f.end);
                w.write(&name, Some(&desc), f.slice(&seq)).unwrap();
            }
        });
    });
    group.finish();
}

/// Per-base mutation throughput at the default and at high error rates.
fn mutate(c: &mut Criterion) {
    let mut group = c.benchmark_group("mutate");
    let seq = random_dna(1_000_000, 6);
    group.throughput(Throughput::Bytes(seq.len() as u64));
    for (name, (s, i, d)) in [
        ("none", (0.0, 0.0, 0.0)),
        ("default", (0.000057, 0.000069, 0.0016)),
        ("high", (0.05, 0.05, 0.05)),
    ] {
        let m = Mutator::new(s, i, d).unwrap();
        group.bench_function(name, |b| {
            let mut rng = StdRng::seed_from_u64(7);
            let mut out = Vec::with_capacity(seq.len() * 2);
            b.iter(|| {
                out.clear();
                m.mutate_into(black_box(&seq), &mut rng, &mut out);
                black_box(&out);
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    sonicate_by_length,
    sonicate_by_depth,
    sonicate_length_model,
    sonicate_to_fasta,
    mutate
);
criterion_main!(benches);
