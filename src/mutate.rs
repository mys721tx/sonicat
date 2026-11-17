//! *In silico* point mutation: substitutions, insertions and deletions.
//!
//! # Model
//!
//! Each input base independently undergoes at most one event:
//!
//! * **substitution** (probability `s`): an `A`, `C`, `G` or `T` is replaced by
//!   one of the other three nucleotides, chosen uniformly; case is preserved.
//!   Other symbols (e.g. `N`) are left unchanged.
//! * **insertion** (probability `i`): the base is kept and a uniformly random
//!   nucleotide is inserted after it.
//! * **deletion** (probability `d`): the base is dropped.
//! * otherwise (probability `1 - s - i - d`) the base is copied unchanged.

use crate::ParamError;
use rand::{Rng, RngExt, distr::weighted::WeightedIndex};

const NUCLEOTIDES: [u8; 4] = *b"ACGT";

/// The fate of a single base under the mutation model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Event {
    Substitution,
    Insertion,
    Deletion,
    Unchanged,
}

/// Events in the order of the weights passed to [`WeightedIndex`].
const EVENTS: [Event; 4] = [
    Event::Substitution,
    Event::Insertion,
    Event::Deletion,
    Event::Unchanged,
];

/// Per-base mutation model. See the [module documentation](self).
#[derive(Debug, Clone)]
pub struct Mutator {
    fate: WeightedIndex<f64>,
}

impl Mutator {
    /// Creates a mutator with per-base substitution, insertion and deletion
    /// probabilities. Each must lie in `[0, 1]` and their sum must not
    /// exceed 1.
    pub fn new(substitution: f64, insertion: f64, deletion: f64) -> Result<Self, ParamError> {
        for (name, p) in [
            ("substitution", substitution),
            ("insertion", insertion),
            ("deletion", deletion),
        ] {
            if !(0.0..=1.0).contains(&p) {
                return Err(ParamError(format!(
                    "{name} probability must be within [0, 1], got {p}"
                )));
            }
        }
        let total = substitution + insertion + deletion;
        // Allow for rounding error in the sum of user-supplied probabilities.
        if total > 1.0 + 1e-12 {
            return Err(ParamError(format!(
                "sum of mutation probabilities must not exceed 1, got {total}"
            )));
        }
        let fate = WeightedIndex::new([substitution, insertion, deletion, (1.0 - total).max(0.0)])
            .map_err(|e| ParamError(e.to_string()))?;
        Ok(Mutator { fate })
    }

    /// Draws the fate of a single base.
    pub(crate) fn sample_event<R: Rng + ?Sized>(&self, rng: &mut R) -> Event {
        EVENTS[rng.sample(&self.fate)]
    }

    /// Mutates a single base, appending zero, one or two bases to `out`.
    pub fn mutate_base<R: Rng + ?Sized>(&self, base: u8, rng: &mut R, out: &mut Vec<u8>) {
        match self.sample_event(rng) {
            Event::Substitution => out.push(substitute(base, rng.random_range(0..3))),
            Event::Insertion => insert(base, rng.random_range(0..4), out),
            Event::Deletion => {}
            Event::Unchanged => out.push(base),
        }
    }

    /// Mutates a whole sequence, appending the result to `out`.
    pub fn mutate_into<R: Rng + ?Sized>(&self, seq: &[u8], rng: &mut R, out: &mut Vec<u8>) {
        out.reserve(seq.len() + seq.len() / 16);
        for &base in seq {
            self.mutate_base(base, rng, out);
        }
    }

    /// Mutates a whole sequence, returning the result.
    pub fn mutate<R: Rng + ?Sized>(&self, seq: &[u8], rng: &mut R) -> Vec<u8> {
        let mut out = Vec::new();
        self.mutate_into(seq, rng, &mut out);
        out
    }
}

/// Replaces a nucleotide with the `pick`-th (`0..3`) of the other three, in
/// `ACGT` order, preserving case. Non-nucleotide symbols are returned
/// unchanged.
pub(crate) fn substitute(base: u8, pick: usize) -> u8 {
    debug_assert!(pick < 3);
    let Some(idx) = NUCLEOTIDES
        .iter()
        .position(|&n| n == base.to_ascii_uppercase())
    else {
        return base;
    };
    // Skip over the original base's index.
    let new = NUCLEOTIDES[if pick >= idx { pick + 1 } else { pick }];
    if base.is_ascii_lowercase() {
        new.to_ascii_lowercase()
    } else {
        new
    }
}

/// Keeps `base` and inserts the `nucleotide`-th (`0..4`) of `ACGT` after it.
pub(crate) fn insert(base: u8, nucleotide: usize, out: &mut Vec<u8>) {
    out.push(base);
    out.push(NUCLEOTIDES[nucleotide]);
}

/// Deterministic tests only: every assertion holds for every possible random
/// draw. Distributional properties are checked in `tests/statistical.rs`.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::{random_dna, rng};
    use std::collections::BTreeSet;

    /// Seeds used to exercise invariants; each must hold for every draw.
    const SEEDS: std::ops::Range<u64> = 0..20;

    // ---- parameter validation ----------------------------------------------

    #[test]
    fn rejects_invalid_parameters() {
        for (s, i, d) in [
            (-0.1, 0.0, 0.0),
            (0.0, -0.1, 0.0),
            (0.0, 0.0, -0.1),
            (1.1, 0.0, 0.0),
            (0.5, 0.5, 0.1),
            (f64::NAN, 0.0, 0.0),
            (0.0, f64::INFINITY, 0.0),
        ] {
            assert!(Mutator::new(s, i, d).is_err(), "accepted s={s} i={i} d={d}");
        }
    }

    #[test]
    fn accepts_boundary_parameters() {
        assert!(Mutator::new(0.0, 0.0, 0.0).is_ok());
        assert!(Mutator::new(1.0, 0.0, 0.0).is_ok());
        assert!(Mutator::new(0.3, 0.3, 0.4).is_ok());
        // Sums that exceed 1 only by floating-point rounding.
        assert!(Mutator::new(0.1, 0.2, 0.7).is_ok());
    }

    // ---- per-draw helpers ------------------------------------------------------

    #[test]
    fn each_rate_drives_its_own_event() {
        // A rate of 1 makes the event certain, so this pins the mapping from
        // parameters to events (the old code had them crossed).
        for (params, event) in [
            ((1.0, 0.0, 0.0), Event::Substitution),
            ((0.0, 1.0, 0.0), Event::Insertion),
            ((0.0, 0.0, 1.0), Event::Deletion),
            ((0.0, 0.0, 0.0), Event::Unchanged),
        ] {
            let m = Mutator::new(params.0, params.1, params.2).unwrap();
            for seed in SEEDS {
                let mut r = rng(seed);
                assert!((0..1_000).all(|_| m.sample_event(&mut r) == event));
            }
        }
    }

    #[test]
    fn substitution_reaches_exactly_the_other_three_nucleotides() {
        for &base in b"ACGTacgt" {
            let targets: BTreeSet<u8> = (0..3).map(|pick| substitute(base, pick)).collect();
            let expected: BTreeSet<u8> = b"ACGT"
                .iter()
                .map(|&n| {
                    if base.is_ascii_lowercase() {
                        n.to_ascii_lowercase()
                    } else {
                        n
                    }
                })
                .filter(|&n| n != base)
                .collect();
            assert_eq!(targets, expected, "base {}", base as char);
        }
    }

    #[test]
    fn substitution_leaves_other_symbols_unchanged() {
        for &base in b"NnRY-*" {
            for pick in 0..3 {
                assert_eq!(substitute(base, pick), base);
            }
        }
    }

    #[test]
    fn insertion_keeps_base_and_can_insert_every_nucleotide() {
        let mut inserted = BTreeSet::new();
        for nucleotide in 0..4 {
            let mut out = Vec::new();
            insert(b'g', nucleotide, &mut out);
            assert_eq!(out.len(), 2);
            assert_eq!(out[0], b'g');
            inserted.insert(out[1]);
        }
        assert_eq!(inserted, b"ACGT".iter().copied().collect());
    }

    // ---- invariants over every draw ------------------------------------------

    #[test]
    fn zero_rates_are_identity() {
        let m = Mutator::new(0.0, 0.0, 0.0).unwrap();
        let seq = random_dna(10_000, 1);
        for seed in SEEDS {
            assert_eq!(m.mutate(&seq, &mut rng(seed)), seq);
            assert_eq!(m.mutate(b"acgtNnRY-*", &mut rng(seed)), b"acgtNnRY-*");
        }
    }

    #[test]
    fn empty_sequence_stays_empty() {
        let m = Mutator::new(0.2, 0.2, 0.2).unwrap();
        for seed in SEEDS {
            assert!(m.mutate(b"", &mut rng(seed)).is_empty());
        }
    }

    #[test]
    fn certain_deletion_removes_everything() {
        let m = Mutator::new(0.0, 0.0, 1.0).unwrap();
        let seq = random_dna(1_000, 4);
        for seed in SEEDS {
            assert!(m.mutate(&seq, &mut rng(seed)).is_empty());
        }
    }

    #[test]
    fn certain_substitution_changes_every_base_and_keeps_case() {
        let m = Mutator::new(1.0, 0.0, 0.0).unwrap();
        let seq = b"ACGTacgtNnRY-";
        for seed in SEEDS {
            let out = m.mutate(seq, &mut rng(seed));
            assert_eq!(out.len(), seq.len());
            for (a, b) in seq[..8].iter().zip(&out[..8]) {
                assert_ne!(a, b);
                assert_eq!(a.is_ascii_lowercase(), b.is_ascii_lowercase());
                assert!(b"ACGT".contains(&b.to_ascii_uppercase()));
            }
            assert_eq!(&out[8..], b"NnRY-");
        }
    }

    #[test]
    fn certain_insertion_doubles_length() {
        let m = Mutator::new(0.0, 1.0, 0.0).unwrap();
        let seq = random_dna(1_000, 9);
        for seed in SEEDS {
            let out = m.mutate(&seq, &mut rng(seed));
            assert_eq!(out.len(), 2 * seq.len());
            let originals: Vec<u8> = out.iter().step_by(2).copied().collect();
            assert_eq!(originals, seq);
            assert!(out.iter().skip(1).step_by(2).all(|b| b"ACGT".contains(b)));
        }
    }

    #[test]
    fn output_length_is_bounded() {
        // Each base yields 0, 1 or 2 bases.
        let m = Mutator::new(0.2, 0.2, 0.2).unwrap();
        let seq = random_dna(1_000, 10);
        for seed in SEEDS {
            assert!(m.mutate(&seq, &mut rng(seed)).len() <= 2 * seq.len());
        }
    }

    #[test]
    fn mutate_into_appends() {
        let m = Mutator::new(0.0, 0.0, 0.0).unwrap();
        let mut out = b"XX".to_vec();
        m.mutate_into(b"ACGT", &mut rng(14), &mut out);
        assert_eq!(out, b"XXACGT");
    }

    #[test]
    fn same_seed_is_reproducible() {
        let m = Mutator::new(0.05, 0.05, 0.05).unwrap();
        let seq = random_dna(10_000, 11);
        for seed in SEEDS {
            assert_eq!(
                m.mutate(&seq, &mut rng(seed)),
                m.mutate(&seq, &mut rng(seed))
            );
        }
    }
}
