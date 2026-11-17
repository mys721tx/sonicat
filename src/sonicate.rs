//! *In silico* sonication: random shearing of a sequence into fragments.
//!
//! # Model
//!
//! For a sequence of length `L`, a sonicator configured with coverage
//! `depth`, mean fragment length `mean` and standard deviation `sd`:
//!
//! 1. draws the number of fragments `N ~ Poisson(depth * L / mean)`, so the
//!    expected number of sequenced bases is `depth * L`;
//! 2. for each fragment, draws a length `ℓ ~ Normal(mean, sd)`, rounded to the
//!    nearest integer and clamped to `[1, L]`;
//! 3. places the fragment at a start position drawn uniformly from
//!    `[0, L - ℓ]`, so every fragment lies entirely within the sequence.
//!
//! An interior base is therefore covered `≈ depth` times on average; bases
//! within one fragment length of either end are covered less, as with real
//! sheared DNA.

use crate::ParamError;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, Normal, Poisson};

/// A fragment `[start, end)` of a source sequence (0-based, half-open).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Fragment {
    pub start: usize,
    pub end: usize,
}

impl Fragment {
    /// Length of the fragment in bases.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Whether the fragment is empty. Fragments produced by a
    /// [`Sonicator`] are never empty.
    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    /// The bases of `seq` covered by this fragment.
    pub fn slice<'a>(&self, seq: &'a [u8]) -> &'a [u8] {
        &seq[self.start..self.end]
    }
}

/// Random fragment generator. See the [module documentation](self) for the
/// model.
#[derive(Debug, Clone)]
pub struct Sonicator {
    depth: f64,
    mean: f64,
    length: Option<Normal<f64>>,
}

impl Sonicator {
    /// Creates a sonicator with per-base coverage `depth`, mean fragment
    /// length `mean` and fragment length standard deviation `sd`.
    ///
    /// `depth` must be finite and non-negative, `mean` finite and at least 1,
    /// and `sd` finite and non-negative. An `sd` of 0 yields fragments of
    /// exactly `mean` (rounded) bases.
    pub fn new(depth: f64, mean: f64, sd: f64) -> Result<Self, ParamError> {
        if !depth.is_finite() || depth < 0.0 {
            return Err(ParamError(format!(
                "depth must be finite and non-negative, got {depth}"
            )));
        }
        if !mean.is_finite() || mean < 1.0 {
            return Err(ParamError(format!(
                "mean fragment length must be finite and at least 1, got {mean}"
            )));
        }
        if !sd.is_finite() || sd < 0.0 {
            return Err(ParamError(format!(
                "fragment length standard deviation must be finite and non-negative, got {sd}"
            )));
        }
        let length = if sd > 0.0 {
            Some(Normal::new(mean, sd).map_err(|e| ParamError(e.to_string()))?)
        } else {
            None
        };
        Ok(Sonicator {
            depth,
            mean,
            length,
        })
    }

    /// Per-base coverage.
    pub fn depth(&self) -> f64 {
        self.depth
    }

    /// Mean fragment length.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Expected number of fragments for a sequence of `seq_len` bases.
    pub fn expected_fragments(&self, seq_len: usize) -> f64 {
        self.depth * seq_len as f64 / self.mean
    }

    /// Draws the length of a single fragment from a sequence of `seq_len`
    /// bases. `seq_len` must be positive.
    fn fragment_length<R: Rng + ?Sized>(&self, seq_len: usize, rng: &mut R) -> usize {
        let raw = match &self.length {
            Some(normal) => normal.sample(rng),
            None => self.mean,
        };
        clamp_length(raw, seq_len)
    }

    /// Shears a sequence of `seq_len` bases into random fragments.
    ///
    /// Fragments are returned in generation order, not sorted by position.
    pub fn fragments<R: Rng + ?Sized>(&self, seq_len: usize, rng: &mut R) -> Vec<Fragment> {
        let lambda = self.expected_fragments(seq_len);
        if seq_len == 0 || lambda <= 0.0 {
            return Vec::new();
        }
        // `lambda` is positive and finite here, and far below
        // `Poisson::MAX_LAMBDA` for any realistic genome.
        let count = Poisson::new(lambda)
            .expect("valid Poisson rate")
            .sample(rng) as usize;

        (0..count)
            .map(|_| {
                let len = self.fragment_length(seq_len, rng);
                let start = rng.random_range(start_range(seq_len, len));
                Fragment {
                    start,
                    end: start + len,
                }
            })
            .collect()
    }
}

/// Converts a raw fragment length draw into a valid length for a sequence of
/// `seq_len` (> 0) bases: rounded to the nearest integer and clamped to
/// `[1, seq_len]`. Negative and NaN draws become 1.
pub(crate) fn clamp_length(raw: f64, seq_len: usize) -> usize {
    // `as` saturates out-of-range floats and maps NaN to 0.
    (raw.round() as usize).clamp(1, seq_len)
}

/// Start positions at which a fragment of `len` bases fits entirely within a
/// sequence of `seq_len` bases. Requires `1 <= len <= seq_len`.
pub(crate) fn start_range(seq_len: usize, len: usize) -> std::ops::RangeInclusive<usize> {
    0..=seq_len - len
}

/// Deterministic tests only: every assertion holds for every possible random
/// draw. Distributional properties are checked in `tests/statistical.rs`.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::{random_dna, rng};

    /// Seeds used to exercise invariants; each must hold for every draw.
    const SEEDS: std::ops::Range<u64> = 0..20;

    // ---- parameter validation ----------------------------------------------

    #[test]
    fn rejects_invalid_parameters() {
        for (depth, mean, sd) in [
            (-1.0, 150.0, 10.0),
            (f64::NAN, 150.0, 10.0),
            (f64::INFINITY, 150.0, 10.0),
            (10.0, 0.0, 10.0),
            (10.0, 0.5, 10.0),
            (10.0, -150.0, 10.0),
            (10.0, f64::NAN, 10.0),
            (10.0, 150.0, -1.0),
            (10.0, 150.0, f64::INFINITY),
        ] {
            assert!(
                Sonicator::new(depth, mean, sd).is_err(),
                "accepted depth={depth} mean={mean} sd={sd}"
            );
        }
    }

    #[test]
    fn accepts_boundary_parameters() {
        assert!(Sonicator::new(0.0, 1.0, 0.0).is_ok());
        assert!(Sonicator::new(1e-9, 1.0, 1e-9).is_ok());
    }

    #[test]
    fn expected_fragments_gives_depth_times_length_over_mean() {
        let s = Sonicator::new(5.0, 100.0, 10.0).unwrap();
        assert_eq!(s.expected_fragments(2_000), 100.0);
        assert_eq!(s.expected_fragments(0), 0.0);
    }

    // ---- per-draw helpers ------------------------------------------------------

    #[test]
    fn clamp_length_rounds_to_nearest() {
        assert_eq!(clamp_length(99.4, 1_000), 99);
        assert_eq!(clamp_length(99.5, 1_000), 100);
        assert_eq!(clamp_length(99.6, 1_000), 100);
    }

    #[test]
    fn clamp_length_never_below_one() {
        for raw in [0.0, 0.4, -0.6, -150.0, f64::NEG_INFINITY, f64::NAN] {
            assert_eq!(clamp_length(raw, 1_000), 1, "raw = {raw}");
        }
    }

    #[test]
    fn clamp_length_never_above_sequence() {
        for raw in [40.6, 150.0, 1e30, f64::INFINITY] {
            assert_eq!(clamp_length(raw, 40), 40, "raw = {raw}");
        }
        assert_eq!(clamp_length(150.0, 1), 1);
    }

    #[test]
    fn start_range_keeps_every_fragment_inside_the_sequence() {
        for seq_len in 1..=60 {
            for len in 1..=seq_len {
                let r = start_range(seq_len, len);
                assert_eq!(*r.start(), 0);
                assert_eq!(r.end() + len, seq_len, "last start ends at the last base");
            }
        }
    }

    // ---- degenerate inputs ---------------------------------------------------

    #[test]
    fn empty_sequence_yields_no_fragments() {
        let s = Sonicator::new(50.0, 150.0, 15.0).unwrap();
        for seed in SEEDS {
            assert!(s.fragments(0, &mut rng(seed)).is_empty());
        }
    }

    #[test]
    fn zero_depth_yields_no_fragments() {
        let s = Sonicator::new(0.0, 150.0, 15.0).unwrap();
        for seed in SEEDS {
            assert!(s.fragments(10_000, &mut rng(seed)).is_empty());
        }
    }

    // ---- invariants over every draw ------------------------------------------

    #[test]
    fn fragments_are_nonempty_and_within_bounds() {
        for &(seq_len, mean, sd) in &[
            (1, 150.0, 30.0),
            (10, 150.0, 15.0),
            (149, 150.0, 15.0),
            (150, 150.0, 15.0),
            (151, 150.0, 15.0),
            (5_000, 150.0, 15.0),
            (5_000, 150.0, 200.0), // sd large enough to draw negatives
            (5_000, 1.0, 0.0),
            (5_000, 4_999.0, 1.0),
        ] {
            let s = Sonicator::new(30.0, mean, sd).unwrap();
            for seed in SEEDS {
                for f in s.fragments(seq_len, &mut rng(seed)) {
                    assert!(!f.is_empty(), "empty fragment {f:?}");
                    assert!(f.end <= seq_len, "{f:?} past end of {seq_len}");
                }
            }
        }
    }

    #[test]
    fn fixed_length_when_sd_is_zero() {
        let s = Sonicator::new(20.0, 150.0, 0.0).unwrap();
        for seed in SEEDS {
            assert!(
                s.fragments(10_000, &mut rng(seed))
                    .iter()
                    .all(|f| f.len() == 150)
            );
        }
    }

    #[test]
    fn sequence_shorter_than_mean_is_returned_whole() {
        let s = Sonicator::new(1_000.0, 150.0, 0.0).unwrap();
        for seed in SEEDS {
            assert!(
                s.fragments(40, &mut rng(seed))
                    .iter()
                    .all(|f| *f == Fragment { start: 0, end: 40 })
            );
        }
    }

    #[test]
    fn fragment_slice_matches_source() {
        let seq = random_dna(2_000, 6);
        let s = Sonicator::new(10.0, 100.0, 20.0).unwrap();
        for seed in SEEDS {
            for f in s.fragments(seq.len(), &mut rng(seed)) {
                assert_eq!(f.slice(&seq), &seq[f.start..f.end]);
                assert_eq!(f.slice(&seq).len(), f.len());
            }
        }
    }

    #[test]
    fn same_seed_is_reproducible() {
        let s = Sonicator::new(30.0, 150.0, 15.0).unwrap();
        for seed in SEEDS {
            assert_eq!(
                s.fragments(10_000, &mut rng(seed)),
                s.fragments(10_000, &mut rng(seed))
            );
        }
    }
}
