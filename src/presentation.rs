use std::fmt;

use crate::word::Word;

/// A power-commutator (PC) presentation of a finite polycyclic group.
///
/// Construction of such presentations happens via [`Builder`](crate::Builder).
///
/// # Definition
/// A polycyclic group has a PC presentation with generators `g0, ..., g{n-1}`,
/// relative orders `p0, ..., p{n-1}`, power relations:
/// - `gi^pi = w_ii` for `0 <= i < n`
///
/// and commutator relations:
/// - `[gi, gj] = gi^-1 * gj^-1 * gi * gj = w_ij` for `0 <= j < i < n`
///
/// where each right-hand side `w_ij` is a normal form [word](Word) in generators `g{j+1}, ..., g{n-1}`.
///
/// # Invariants
/// - `relative_orders.len() == n` and each `p_i >= 2` is prime.
/// - `powers.len() == n`.
/// - `commutators.len() == num_commutators(num_gens) = (n * (n - 1)) / 2`.
/// - Commutators are stored in a flat lower-triangular sequence indexed by [`commutator_index`](Self::commutator_index).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Presentation {
    /// The relative orders (prime orders of the cyclic factor groups) `p_0, ..., p_{n-1}`.
    relative_orders: Vec<u32>,

    /// Power relations `g_i^{p_i}` for each `i = 0, ..., n-1`.
    powers: Vec<Word>,

    /// Commutator relations `[g_i, g_j]` for `0 <= j < i < n`.
    /// Stored as a flattened lower-triangular matrix for optimal cache locality.
    commutators: Vec<Word>,
}

impl Presentation {
    /// Computes the number of commutator pairs `(i, j)` with `0 <= j < i < n`, which is `n * (n - 1) / 2`.
    #[inline]
    pub const fn num_commutators(num_gens: usize) -> usize {
        (num_gens * num_gens.saturating_sub(1)) / 2
    }

    /// Computes the flat vector index for the commutator `[g_i, g_j]` where `i > j`.
    ///
    /// The indexing follows standard row-major lower-triangular layout:
    /// `index(i, j) = (i * (i - 1)) / 2 + j`.
    ///
    /// # Panics
    /// Panics in debug mode if `i <= j`.
    #[inline]
    pub const fn commutator_index(i: usize, j: usize) -> usize {
        debug_assert!(i > j, "Commutator indexing requires i > j");
        (i * (i - 1)) / 2 + j
    }

    /// Creates a new `Presentation` from its constituent parts.
    ///
    /// Expected to be called by [`Builder`](crate::Builder).
    #[inline]
    pub(crate) fn new(
        relative_orders: Vec<u32>,
        powers: Vec<Word>,
        commutators: Vec<Word>,
    ) -> Self {
        debug_assert_eq!(powers.len(), relative_orders.len());
        debug_assert_eq!(
            commutators.len(),
            Self::num_commutators(relative_orders.len())
        );
        Self {
            relative_orders,
            powers,
            commutators,
        }
    }

    /// Returns the total number of generators `n` in the presentation.
    #[inline]
    pub fn num_gens(&self) -> usize {
        self.relative_orders.len()
    }

    /// Returns the relative order `p_i` for generator `g_i`.
    ///
    /// # Panics
    /// Panics if `i >= self.num_gens()`.
    #[inline]
    pub fn relative_order(&self, i: usize) -> u32 {
        self.relative_orders[i]
    }

    /// Returns a slice of all generator relative orders `[p_0, ..., p_{n-1}]`.
    #[inline]
    pub fn relative_orders(&self) -> &[u32] {
        &self.relative_orders
    }

    /// Returns a reference to the normal form [word](Word) for the power relation `g_i^{p_i}`.
    ///
    /// # Panics
    /// Panics if `i >= self.num_gens()`.
    #[inline]
    pub fn power(&self, i: usize) -> &Word {
        &self.powers[i]
    }

    /// Returns a slice of all power relation words `[g_0^{p_0}, ..., g_{n-1}^{p_{n-1}}]`.
    #[inline]
    pub fn powers(&self) -> &[Word] {
        &self.powers
    }

    /// Returns a reference to the normal form [word](Word) for the commutator relation `[g_i, g_j]`.
    ///
    /// # Panics
    /// Panics if `i <= j` or if `i >= self.num_gens()`.
    #[inline]
    pub fn commutator(&self, i: usize, j: usize) -> &Word {
        assert!(
            i > j,
            "Commutator relation requires i > j, but got i = {i}, j = {j}"
        );
        &self.commutators[Self::commutator_index(i, j)]
    }

    /// Returns a slice of all commutator relation words.
    #[inline]
    pub fn commutators(&self) -> &[Word] {
        &self.commutators
    }

    /// Computes the total group order `p0 * p1 * ... * p_{n-1}`.
    #[inline]
    pub fn order(&self) -> u128 {
        self.relative_orders
            .iter()
            .fold(1u128, |acc, &p| acc * p as u128)
    }
}

impl fmt::Display for Presentation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.num_gens();
        if n == 0 {
            return write!(f, "< >");
        }

        write!(f, "< ")?;

        // Format generators: g0, g1, ..., g{n-1}
        for i in 0..n {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "g{i}")?;
        }

        write!(f, " | ")?;

        let mut first_rel = true;

        // Power relations: gi^pi = word
        for i in 0..n {
            if !first_rel {
                write!(f, ", ")?;
            }
            first_rel = false;
            let p = self.relative_order(i);
            let power_word = self.power(i);
            write!(f, "g{i}^{p} = {power_word}")?;
        }

        // Commutator relations: [gi, gj] = word
        for i in 1..n {
            for j in 0..i {
                if !first_rel {
                    write!(f, ", ")?;
                }
                first_rel = false;
                let comm_word = self.commutator(i, j);
                write!(f, "[g{i}, g{j}] = {comm_word}")?;
            }
        }

        write!(f, " >")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_presentation_display() {
        let empty_pres = Presentation::new(Vec::new(), Vec::new(), Vec::new());
        assert_eq!(empty_pres.to_string(), "< >");

        let rel_orders = vec![2, 3];
        let powers = vec![Word::identity(), Word::from_term(0, 1)];
        let commutators = vec![Word::from_term(1, 2)]; // [g1, g0]

        let pres = Presentation::new(rel_orders, powers, commutators);
        assert_eq!(
            pres.to_string(),
            "< g0, g1 | g0^2 = 1, g1^3 = g0, [g1, g0] = g1^2 >"
        );
    }

    #[test]
    fn test_num_commutators() {
        assert_eq!(Presentation::num_commutators(0), 0);
        assert_eq!(Presentation::num_commutators(1), 0);
        assert_eq!(Presentation::num_commutators(2), 1);
        assert_eq!(Presentation::num_commutators(3), 3);
        assert_eq!(Presentation::num_commutators(4), 6);
        assert_eq!(Presentation::num_commutators(5), 10);
    }

    #[test]
    fn test_commutator_indexing() {
        // Lower-triangular indexing: (i, j) with i > j
        assert_eq!(Presentation::commutator_index(1, 0), 0);
        assert_eq!(Presentation::commutator_index(2, 0), 1);
        assert_eq!(Presentation::commutator_index(2, 1), 2);
        assert_eq!(Presentation::commutator_index(3, 0), 3);
        assert_eq!(Presentation::commutator_index(3, 1), 4);
        assert_eq!(Presentation::commutator_index(3, 2), 5);
    }

    #[test]
    fn test_empty_presentation() {
        let pres = Presentation::new(Vec::new(), Vec::new(), Vec::new());
        assert_eq!(pres.num_gens(), 0);
        assert_eq!(pres.order(), 1);
        assert!(pres.relative_orders().is_empty());
        assert!(pres.powers().is_empty());
        assert!(pres.commutators().is_empty());
    }

    #[test]
    fn test_presentation_construction_and_getters() {
        let rel_orders = vec![2, 3, 2];
        let powers = vec![Word::identity(), Word::from_term(2, 1), Word::identity()];
        let commutators = vec![
            Word::from_term(2, 1), // [g1, g0] -> index 0
            Word::identity(),      // [g2, g0] -> index 1
            Word::identity(),      // [g2, g1] -> index 2
        ];

        let pres = Presentation::new(rel_orders.clone(), powers.clone(), commutators.clone());

        assert_eq!(pres.num_gens(), 3);
        assert_eq!(pres.order(), 2 * 3 * 2);
        assert_eq!(pres.relative_orders(), &[2, 3, 2]);
        assert_eq!(pres.relative_order(0), 2);
        assert_eq!(pres.relative_order(1), 3);
        assert_eq!(pres.relative_order(2), 2);

        assert_eq!(pres.powers(), &powers[..]);
        assert_eq!(pres.power(1), &Word::from_term(2, 1));

        assert_eq!(pres.commutators(), &commutators[..]);
        assert_eq!(pres.commutator(1, 0), &Word::from_term(2, 1));
        assert_eq!(pres.commutator(2, 0), &Word::identity());
        assert_eq!(pres.commutator(2, 1), &Word::identity());

        // Test Clone and PartialEq
        let pres_clone = pres.clone();
        assert_eq!(pres, pres_clone);
    }

    #[test]
    #[should_panic(expected = "Commutator relation requires i > j")]
    fn test_invalid_commutator_access() {
        let pres = Presentation::new(
            vec![2, 2],
            vec![Word::identity(); 2],
            vec![Word::identity()],
        );
        let _ = pres.commutator(0, 1);
    }
}
