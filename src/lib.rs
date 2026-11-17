//! Core library for **Sonicat**, a suite of tools for *in silico* sonication
//! and sequence manipulation.
//!
//! All randomness is injected through a caller-supplied [`rand::Rng`], so the
//! simulators are reproducible when driven by a seeded generator.

pub mod mutate;
pub mod sonicate;

pub use mutate::Mutator;
pub use sonicate::{Fragment, Sonicator};

use std::fmt;

/// Error raised when a simulator is configured with invalid parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct ParamError(pub String);

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ParamError {}

#[cfg(test)]
pub(crate) mod test_util {
    use rand::{RngExt, SeedableRng, rngs::StdRng};

    pub fn rng(seed: u64) -> StdRng {
        StdRng::seed_from_u64(seed)
    }

    /// A random ACGT sequence of length `n`.
    pub fn random_dna(n: usize, seed: u64) -> Vec<u8> {
        let mut rng = rng(seed);
        (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
    }
}
