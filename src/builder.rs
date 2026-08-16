use std::fmt;

use crate::util::is_prime;
use crate::{Presentation, Word};

/// Errors that can occur during the construction of a [presentation](Presentation) via [`Builder`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuilderError {
    /// The specified relative order is not a prime number (`order >= 2`).
    InvalidRelativeOrder { order: u32 },

    /// The specified generator index does not exist in this presentation.
    GeneratorOutOfBounds { index: usize, max: usize },

    /// Violates the PC condition: relations defining `gi` or `[gi, gj]` must only contain generators `gk` where `k > i` (for powers) or `k > j` (for commutators).
    InvalidPcCondition {
        relation_for: usize,
        invalid_generator: usize,
    },

    /// Attempted to define a commutator `[gi, gj]` where `i <= j`.
    InvalidCommutatorPair { i: usize, j: usize },
}

impl fmt::Display for BuilderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRelativeOrder { order } => {
                write!(
                    f,
                    "Invalid relative order {order}: must be a prime number >= 2"
                )
            }
            Self::GeneratorOutOfBounds { index, max } => {
                write!(
                    f,
                    "Generator index {index} is out of bounds (maximum valid index is {max})"
                )
            }
            Self::InvalidPcCondition {
                relation_for,
                invalid_generator,
            } => {
                write!(
                    f,
                    "PC condition violated: relation for g{relation_for} contains generator g{invalid_generator} (must be strictly greater)"
                )
            }
            Self::InvalidCommutatorPair { i, j } => {
                write!(f, "Invalid commutator pair [g{i}, g{j}]: requires i > j")
            }
        }
    }
}

impl std::error::Error for BuilderError {}

/// A builder to ergonomically construct a [PC presentation](Presentation).
///
/// Unspecified power and commutator relations are implicitly trivial (the identity [word](Word)).
///
/// # Example: Constructing S_3
/// ```
/// use pcgroup::{Builder, Word};
///
/// // Construct the symmetric group S_3 (dihedral group D_3 of order 6):
/// // < g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >
/// let pres = Builder::new(vec![2, 3])
///     .unwrap()
///     .add_commutator(1, 0, Word::from_term(1, 1))
///     .unwrap()
///     .build();
///
/// assert_eq!(pres.num_gens(), 2);
/// assert_eq!(pres.order(), 6);
/// assert_eq!(
///     pres.to_string(),
///     "< g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >"
/// );
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Builder {
    relative_orders: Vec<u32>,
    powers: Vec<Word>,
    commutators: Vec<Word>,
}

impl Builder {
    /// Starts building a presentation with the given prime relative orders.
    ///
    /// # Errors
    /// Returns [`BuilderError::InvalidRelativeOrder`] if any given order is not a prime number (`>= 2`).
    pub fn new(relative_orders: Vec<u32>) -> Result<Self, BuilderError> {
        for &order in &relative_orders {
            if !is_prime(order) {
                return Err(BuilderError::InvalidRelativeOrder { order });
            }
        }

        let num_gens = relative_orders.len();
        let num_comms = Presentation::num_commutators(num_gens);

        Ok(Self {
            relative_orders,
            powers: vec![Word::identity(); num_gens],
            commutators: vec![Word::identity(); num_comms],
        })
    }

    /// Returns the number of generators in the presentation being built.
    #[inline]
    pub fn num_gens(&self) -> usize {
        self.relative_orders.len()
    }

    /// Validates that a word only contains valid generators according to the PC condition.
    fn validate_word(&self, relation_for: usize, word: &Word) -> Result<(), BuilderError> {
        let num_gens = self.num_gens();
        let max = num_gens.saturating_sub(1);
        for term in word.iter() {
            if term.gen_index >= num_gens {
                return Err(BuilderError::GeneratorOutOfBounds {
                    index: term.gen_index,
                    max,
                });
            }
            if term.gen_index <= relation_for {
                return Err(BuilderError::InvalidPcCondition {
                    relation_for,
                    invalid_generator: term.gen_index,
                });
            }
        }
        Ok(())
    }

    /// Adds or overwrites a power relation: `gi^pi = word`.
    ///
    /// # Errors
    /// Returns [`BuilderError::GeneratorOutOfBounds`] if `i >= num_gens` or if `word` contains out-of-bounds generators.
    /// Returns [`BuilderError::InvalidPcCondition`] if `word` contains generators `<= i`.
    pub fn add_power(mut self, i: usize, word: Word) -> Result<Self, BuilderError> {
        let num_gens = self.num_gens();
        if i >= num_gens {
            return Err(BuilderError::GeneratorOutOfBounds {
                index: i,
                max: num_gens.saturating_sub(1),
            });
        }
        self.validate_word(i, &word)?;

        self.powers[i] = word;
        Ok(self)
    }

    /// Adds or overwrites a commutator relation: `[gi, gj] = word` for `i > j`.
    ///
    /// # Errors
    /// Returns [`BuilderError::GeneratorOutOfBounds`] if `i >= num_gens`, `j >= num_gens`, or if `word` contains out-of-bounds generators.
    /// Returns [`BuilderError::InvalidCommutatorPair`] if `i <= j`.
    /// Returns [`BuilderError::InvalidPcCondition`] if `word` contains generators `<= j`.
    pub fn add_commutator(mut self, i: usize, j: usize, word: Word) -> Result<Self, BuilderError> {
        let num_gens = self.num_gens();
        let max = num_gens.saturating_sub(1);
        if i >= num_gens {
            return Err(BuilderError::GeneratorOutOfBounds { index: i, max });
        }
        if j >= num_gens {
            return Err(BuilderError::GeneratorOutOfBounds { index: j, max });
        }
        if i <= j {
            return Err(BuilderError::InvalidCommutatorPair { i, j });
        }

        // The commutator [gi, gj] moves gj past gi, landing in G_{j+1} (generators > j)
        self.validate_word(j, &word)?;

        let flat_idx = Presentation::commutator_index(i, j);
        self.commutators[flat_idx] = word;
        Ok(self)
    }

    /// Finalizes the presentation, returning an immutable [`Presentation`].
    pub fn build(self) -> Presentation {
        Presentation::new(self.relative_orders, self.powers, self.commutators)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builder_s3() {
        let pres = Builder::new(vec![2, 3])
            .unwrap()
            .add_commutator(1, 0, Word::from_term(1, 1))
            .unwrap()
            .build();

        assert_eq!(pres.num_gens(), 2);
        assert_eq!(pres.order(), 6);
        assert_eq!(
            pres.to_string(),
            "< g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >"
        );
    }

    #[test]
    fn test_builder_klein4() {
        let pres = Builder::new(vec![2, 2]).unwrap().build();
        assert_eq!(pres.num_gens(), 2);
        assert_eq!(pres.order(), 4);
        assert_eq!(
            pres.to_string(),
            "< g0, g1 | g0^2 = 1, g1^2 = 1, [g1, g0] = 1 >"
        );
    }

    #[test]
    fn test_builder_c4() {
        let pres = Builder::new(vec![2, 2])
            .unwrap()
            .add_power(0, Word::from_term(1, 1))
            .unwrap()
            .build();

        assert_eq!(pres.num_gens(), 2);
        assert_eq!(pres.order(), 4);
        assert_eq!(
            pres.to_string(),
            "< g0, g1 | g0^2 = g1, g1^2 = 1, [g1, g0] = 1 >"
        );
    }

    #[test]
    fn test_invalid_relative_orders() {
        assert_eq!(
            Builder::new(vec![0]),
            Err(BuilderError::InvalidRelativeOrder { order: 0 })
        );
        assert_eq!(
            Builder::new(vec![1]),
            Err(BuilderError::InvalidRelativeOrder { order: 1 })
        );
        assert_eq!(
            Builder::new(vec![4]),
            Err(BuilderError::InvalidRelativeOrder { order: 4 })
        );
        assert_eq!(
            Builder::new(vec![2, 9]),
            Err(BuilderError::InvalidRelativeOrder { order: 9 })
        );
    }

    #[test]
    fn test_generator_out_of_bounds() {
        let builder = Builder::new(vec![2, 3]).unwrap();

        // Out of bounds power generator index
        assert_eq!(
            builder.clone().add_power(2, Word::identity()),
            Err(BuilderError::GeneratorOutOfBounds { index: 2, max: 1 })
        );

        // Out of bounds word generator inside power
        assert_eq!(
            builder.clone().add_power(0, Word::from_term(3, 1)),
            Err(BuilderError::GeneratorOutOfBounds { index: 3, max: 1 })
        );

        // Out of bounds commutator index
        assert_eq!(
            builder.clone().add_commutator(2, 0, Word::identity()),
            Err(BuilderError::GeneratorOutOfBounds { index: 2, max: 1 })
        );
    }

    #[test]
    fn test_invalid_pc_condition() {
        let builder = Builder::new(vec![2, 3, 2]).unwrap();

        // Power relation for g0 containing g0 (must be > 0)
        assert_eq!(
            builder.clone().add_power(0, Word::from_term(0, 1)),
            Err(BuilderError::InvalidPcCondition {
                relation_for: 0,
                invalid_generator: 0,
            })
        );

        // Commutator [g2, g1] containing g1 (must be > 1)
        assert_eq!(
            builder.clone().add_commutator(2, 1, Word::from_term(1, 1)),
            Err(BuilderError::InvalidPcCondition {
                relation_for: 1,
                invalid_generator: 1,
            })
        );
    }

    #[test]
    fn test_invalid_commutator_pair() {
        let builder = Builder::new(vec![2, 3]).unwrap();

        assert_eq!(
            builder.clone().add_commutator(0, 1, Word::identity()),
            Err(BuilderError::InvalidCommutatorPair { i: 0, j: 1 })
        );
        assert_eq!(
            builder.add_commutator(1, 1, Word::identity()),
            Err(BuilderError::InvalidCommutatorPair { i: 1, j: 1 })
        );
    }

    #[test]
    fn test_error_display() {
        let err1 = BuilderError::InvalidRelativeOrder { order: 6 };
        assert_eq!(
            err1.to_string(),
            "Invalid relative order 6: must be a prime number >= 2"
        );

        let err2 = BuilderError::GeneratorOutOfBounds { index: 3, max: 2 };
        assert_eq!(
            err2.to_string(),
            "Generator index 3 is out of bounds (maximum valid index is 2)"
        );

        let err3 = BuilderError::InvalidPcCondition {
            relation_for: 1,
            invalid_generator: 1,
        };
        assert_eq!(
            err3.to_string(),
            "PC condition violated: relation for g1 contains generator g1 (must be strictly greater)"
        );

        let err4 = BuilderError::InvalidCommutatorPair { i: 0, j: 1 };
        assert_eq!(
            err4.to_string(),
            "Invalid commutator pair [g0, g1]: requires i > j"
        );
    }
}
