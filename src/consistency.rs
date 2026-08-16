use std::fmt;

use crate::{Collector, Presentation};

/// Errors that can occur when verifying the consistency of a [PC presentation](Presentation).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConsistencyError {
    /// Inconsistency in the overlap of three distinct generators `gk * gj * gi` where `k > j > i`.
    ThreeGenerators {
        /// The largest generator index `k`.
        k: usize,
        /// The middle generator index `j`.
        j: usize,
        /// The smallest generator index `i`.
        i: usize,
    },

    /// Inconsistency in the overlap of a power and a smaller generator `gj^pj * gi` where `j > i`.
    PowerAndSmallerGenerator {
        /// The larger generator index `j`.
        j: usize,
        /// The smaller generator index `i`.
        i: usize,
    },

    /// Inconsistency in the overlap of a larger generator and a power `gk * gj^pj` where `k > j`.
    LargerGeneratorAndPower {
        /// The larger generator index `k`.
        k: usize,
        /// The smaller generator index `j`.
        j: usize,
    },

    /// Inconsistency in the overlap of a power and itself `gi^(pi + 1)` for generator `gi`.
    PowerAndItself {
        /// The generator index `i`.
        i: usize,
    },
}

impl fmt::Display for ConsistencyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ThreeGenerators { k, j, i } => {
                write!(
                    f,
                    "Consistency check failed for three generators overlap: (g{k} * g{j}) * g{i} != g{k} * (g{j} * g{i})"
                )
            }
            Self::PowerAndSmallerGenerator { j, i } => {
                write!(
                    f,
                    "Consistency check failed for power and smaller generator overlap: (g{j}^p{j}) * g{i} != g{j}^(p{j}-1) * (g{j} * g{i})"
                )
            }
            Self::LargerGeneratorAndPower { k, j } => {
                write!(
                    f,
                    "Consistency check failed for larger generator and power overlap: (g{k} * g{j}) * g{j}^(p{j}-1) != g{k} * (g{j}^p{j})"
                )
            }
            Self::PowerAndItself { i } => {
                write!(
                    f,
                    "Consistency check failed for power and itself overlap: (g{i}^p{i}) * g{i} != g{i} * (g{i}^p{i})"
                )
            }
        }
    }
}

impl std::error::Error for ConsistencyError {}

/// Verifies the consistency (confluence) of a polycyclic (PC) presentation in `O(n^3)` time,
/// where `n` is the number of generators.
///
/// It is assumed that the basic structural validity conditions checked by [`Builder`](crate::Builder)
/// are satisfied (i.e. prime relative orders and PC condition on relation words).
///
/// A PC presentation is consistent if and only if normal forms are unique, which occurs if all
/// critical word overlaps evaluate to identical normal forms under both left-associated and
/// right-associated rewritings:
///
/// 1. **Three distinct generators** (`k > j > i`):
///    - Left: `(gk * gj) * gi -> (gj * gk * [gk, gj]) * gi`
///    - Right: `gk * (gj * gi) -> gk * (gi * gj * [gj, gi])`
/// 2. **Power and smaller generator** (`j > i`):
///    - Left: `(gj^pj) * gi -> w_jj * gi`
///    - Right: `gj^(pj - 1) * (gj * gi) -> gj^(pj - 1) * (gi * gj * [gj, gi])`
/// 3. **Larger generator and power** (`k > j`):
///    - Left: `(gk * gj) * gj^(pj - 1) -> (gj * gk * [gk, gj]) * gj^(pj - 1)`
///    - Right: `gk * (gj^pj) -> gk * w_jj`
/// 4. **Power and itself** (`0 <= i < n`):
///    - Left: `(gi^pi) * gi -> w_ii * gi`
///    - Right: `gi * (gi^pi) -> gi * w_ii`
///
/// # Errors
/// Returns a [`ConsistencyError`] describing the first critical overlap that fails to collect
/// to the same normal form.
pub fn verify_consistency(pres: &Presentation) -> Result<(), ConsistencyError> {
    let n = pres.num_gens();

    // 1. Overlap of three distinct generators (k > j > i)
    for k in 2..n {
        for j in 1..k {
            for i in 0..j {
                let mut col_left = Collector::new(pres);
                col_left.collect_generator(j, 1);
                col_left.collect_generator(k, 1);
                col_left.collect(pres.commutator(k, j));
                col_left.collect_generator(i, 1);

                let mut col_right = Collector::new(pres);
                col_right.collect_generator(k, 1);
                col_right.collect_generator(i, 1);
                col_right.collect_generator(j, 1);
                col_right.collect(pres.commutator(j, i));

                if col_left.into_word() != col_right.into_word() {
                    return Err(ConsistencyError::ThreeGenerators { k, j, i });
                }
            }
        }
    }

    // 2. Overlap of a power and a smaller generator (j > i)
    for j in 1..n {
        let p_j = pres.relative_order(j);
        for i in 0..j {
            let mut col_left = Collector::new(pres);
            col_left.collect(pres.power(j));
            col_left.collect_generator(i, 1);

            let mut col_right = Collector::new(pres);
            col_right.collect_generator(j, p_j - 1);
            col_right.collect_generator(i, 1);
            col_right.collect_generator(j, 1);
            col_right.collect(pres.commutator(j, i));

            if col_left.into_word() != col_right.into_word() {
                return Err(ConsistencyError::PowerAndSmallerGenerator { j, i });
            }
        }
    }

    // 3. Overlap of a larger generator and a power (k > j)
    for k in 1..n {
        for j in 0..k {
            let p_j = pres.relative_order(j);
            let mut col_left = Collector::new(pres);
            col_left.collect_generator(j, 1);
            col_left.collect_generator(k, 1);
            col_left.collect(pres.commutator(k, j));
            col_left.collect_generator(j, p_j - 1);

            let mut col_right = Collector::new(pres);
            col_right.collect_generator(k, 1);
            col_right.collect(pres.power(j));

            if col_left.into_word() != col_right.into_word() {
                return Err(ConsistencyError::LargerGeneratorAndPower { k, j });
            }
        }
    }

    // 4. Overlap of a power and itself (for all i)
    for i in 0..n {
        let mut col_left = Collector::new(pres);
        col_left.collect(pres.power(i));
        col_left.collect_generator(i, 1);

        let mut col_right = Collector::new(pres);
        col_right.collect_generator(i, 1);
        col_right.collect(pres.power(i));

        if col_left.into_word() != col_right.into_word() {
            return Err(ConsistencyError::PowerAndItself { i });
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Builder, Word};

    #[test]
    fn test_error_three_generators() {
        // 4 generators of order 2
        // Relations: [g2, g1] = g3, [g3, g0] = g3
        let pres = Builder::new(vec![2, 2, 2, 2])
            .unwrap()
            .add_commutator(2, 1, Word::from_term(3, 1))
            .unwrap()
            .add_commutator(3, 0, Word::from_term(3, 1))
            .unwrap()
            .build();

        assert_eq!(
            verify_consistency(&pres),
            Err(ConsistencyError::ThreeGenerators { k: 2, j: 1, i: 0 })
        );
    }

    #[test]
    fn test_error_power_and_smaller_generator() {
        // 3 generators of order 2
        // Relations: g1^2 = g2, [g2, g0] = g2
        let pres = Builder::new(vec![2, 2, 2])
            .unwrap()
            .add_power(1, Word::from_term(2, 1))
            .unwrap()
            .add_commutator(2, 0, Word::from_term(2, 1))
            .unwrap()
            .build();

        assert_eq!(
            verify_consistency(&pres),
            Err(ConsistencyError::PowerAndSmallerGenerator { j: 1, i: 0 })
        );
    }

    #[test]
    fn test_error_larger_generator_and_power() {
        // 2 generators of order 2
        // Relations: g0^2 = g1, [g1, g0] = g1
        let pres = Builder::new(vec![2, 2])
            .unwrap()
            .add_power(0, Word::from_term(1, 1))
            .unwrap()
            .add_commutator(1, 0, Word::from_term(1, 1))
            .unwrap()
            .build();

        assert_eq!(
            verify_consistency(&pres),
            Err(ConsistencyError::LargerGeneratorAndPower { k: 1, j: 0 })
        );
    }

    #[test]
    fn test_error_power_and_itself() {
        // 2 generators with orders [2, 3]
        // Relations: g0^2 = g1, [g1, g0] = g1
        let pres = Builder::new(vec![2, 3])
            .unwrap()
            .add_power(0, Word::from_term(1, 1))
            .unwrap()
            .add_commutator(1, 0, Word::from_term(1, 1))
            .unwrap()
            .build();

        assert_eq!(
            verify_consistency(&pres),
            Err(ConsistencyError::PowerAndItself { i: 0 })
        );
    }

    #[test]
    fn test_consistent_trivial_group() {
        let pres = Builder::new(vec![]).unwrap().build();
        assert_eq!(verify_consistency(&pres), Ok(()));
    }

    #[test]
    fn test_consistent_cyclic_group_c3() {
        let pres = Builder::new(vec![3]).unwrap().build();
        assert_eq!(verify_consistency(&pres), Ok(()));
    }

    #[test]
    fn test_consistent_quaternion_group_q8() {
        // Q_8: < g0, g1, g2 | g0^2 = g2, g1^2 = g2, g2^2 = 1, [g1, g0] = g2, [g2, g0] = 1, [g2, g1] = 1 >
        let pres = Builder::new(vec![2, 2, 2])
            .unwrap()
            .add_power(0, Word::from_term(2, 1))
            .unwrap()
            .add_power(1, Word::from_term(2, 1))
            .unwrap()
            .add_commutator(1, 0, Word::from_term(2, 1))
            .unwrap()
            .build();

        assert_eq!(verify_consistency(&pres), Ok(()));
    }
}
