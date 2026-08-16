use std::fmt;

/// Represents a single generator-exponent term `gi^e` of a normal form [word](Word) implicitly associated to some [PC presentation](Presentation).
/// The following invariants must hold (but are not enforced by the type):
///
/// # Invariants
/// - `i = gen_index` is a valid 0-based generator index (`0 <= i < n`) where `n` is the [number of generators](Presentation::num_gens).
/// - `e = exponent` is strictly positive (`e > 0`) and bounded by `p_i` where `p_i` is the [relative order](Presentation::relative_order) of the `i`-th generator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Term {
    /// The generator index `i`.
    pub gen_index: usize,

    /// The exponent `e`.
    pub exponent: u32,
}

impl Term {
    /// Creates a new term with the given generator index and exponent.
    #[inline]
    pub const fn new(gen_index: usize, exponent: u32) -> Self {
        Self {
            gen_index,
            exponent,
        }
    }
}

impl From<(usize, u32)> for Term {
    #[inline]
    fn from((gen_index, exponent): (usize, u32)) -> Self {
        Self::new(gen_index, exponent)
    }
}

impl From<Term> for (usize, u32) {
    #[inline]
    fn from(term: Term) -> Self {
        (term.gen_index, term.exponent)
    }
}

impl fmt::Display for Term {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.exponent == 1 {
            write!(f, "g{}", self.gen_index)
        } else {
            write!(f, "g{}^{}", self.gen_index, self.exponent)
        }
    }
}

/// Represents a word in normal form associated with a [PC presentation](Presentation).
/// An empty sequence of [terms](Term) represents the group identity (`1`).
///
/// # Invariants
/// On top of the invariants of [Term] the sequence of generator indices must be strictly increasing: `i1 < i2 < ... < ik`.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Word {
    /// A vector of terms with strictly increasing `gen_index`.
    pub terms: Vec<Term>,
}

impl Word {
    /// Creates an empty normal form word representing the group identity.
    #[inline]
    pub const fn identity() -> Self {
        Self { terms: Vec::new() }
    }

    /// Creates a new `Word` from a vector of terms.
    #[inline]
    pub const fn new(terms: Vec<Term>) -> Self {
        Self { terms }
    }

    /// Creates a `Word` consisting of a single term $g_i^e$.
    ///
    /// If `exponent == 0`, returns the identity word.
    #[inline]
    pub fn from_term(gen_index: usize, exponent: u32) -> Self {
        if exponent == 0 {
            Self::identity()
        } else {
            Self {
                terms: vec![Term::new(gen_index, exponent)],
            }
        }
    }

    /// Returns `true` if this word represents the group identity.
    #[inline]
    pub fn is_identity(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns `true` if the word has no terms (equivalent to [`is_identity`](Self::is_identity)).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the number of terms in this word.
    #[inline]
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Returns a slice over the terms of this word.
    #[inline]
    pub fn terms(&self) -> &[Term] {
        &self.terms
    }

    /// Returns an iterator over references to the terms of this word.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, Term> {
        self.terms.iter()
    }

    /// Returns the leading term if non-identity.
    #[inline]
    pub fn leading_term(&self) -> Option<Term> {
        self.terms.first().copied()
    }

    /// Returns the leading generator index if non-identity.
    #[inline]
    pub fn leading_generator(&self) -> Option<usize> {
        self.terms.first().map(|term| term.gen_index)
    }

    /// Returns the leading generator exponent if non-identity.
    #[inline]
    pub fn leading_exponent(&self) -> Option<u32> {
        self.terms.first().map(|term| term.exponent)
    }
}

impl fmt::Display for Word {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_identity() {
            write!(f, "1")
        } else {
            for (i, term) in self.terms.iter().enumerate() {
                if i > 0 {
                    write!(f, " ")?;
                }
                write!(f, "{}", term)?;
            }
            Ok(())
        }
    }
}

impl From<Vec<Term>> for Word {
    #[inline]
    fn from(terms: Vec<Term>) -> Self {
        Self { terms }
    }
}

impl From<&[Term]> for Word {
    #[inline]
    fn from(terms: &[Term]) -> Self {
        Self {
            terms: terms.to_vec(),
        }
    }
}

impl FromIterator<Term> for Word {
    #[inline]
    fn from_iter<I: IntoIterator<Item = Term>>(iter: I) -> Self {
        Self {
            terms: iter.into_iter().collect(),
        }
    }
}

impl IntoIterator for Word {
    type Item = Term;
    type IntoIter = std::vec::IntoIter<Term>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.terms.into_iter()
    }
}

impl<'a> IntoIterator for &'a Word {
    type Item = &'a Term;
    type IntoIter = std::slice::Iter<'a, Term>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.terms.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_term_creation_and_display() {
        let t1 = Term::new(0, 1);
        assert_eq!(t1.gen_index, 0);
        assert_eq!(t1.exponent, 1);
        assert_eq!(t1.to_string(), "g0");

        let t2 = Term::new(3, 4);
        assert_eq!(t2.gen_index, 3);
        assert_eq!(t2.exponent, 4);
        assert_eq!(t2.to_string(), "g3^4");
    }

    #[test]
    fn test_term_conversions() {
        let t = Term::from((2, 5));
        assert_eq!(t, Term::new(2, 5));
        let tuple: (usize, u32) = t.into();
        assert_eq!(tuple, (2, 5));
    }

    #[test]
    fn test_word_identity_and_default() {
        let w_id = Word::identity();
        let w_def = Word::default();
        assert_eq!(w_id, w_def);
        assert!(w_id.is_identity());
        assert!(w_id.is_empty());
        assert_eq!(w_id.len(), 0);
        assert_eq!(w_id.leading_term(), None);
        assert_eq!(w_id.leading_generator(), None);
        assert_eq!(w_id.leading_exponent(), None);
        assert_eq!(w_id.to_string(), "1");
    }

    #[test]
    fn test_word_from_term() {
        let w1 = Word::from_term(1, 3);
        assert!(!w1.is_identity());
        assert_eq!(w1.len(), 1);
        assert_eq!(w1.leading_term(), Some(Term::new(1, 3)));
        assert_eq!(w1.leading_generator(), Some(1));
        assert_eq!(w1.leading_exponent(), Some(3));
        assert_eq!(w1.to_string(), "g1^3");

        let w_zero = Word::from_term(1, 0);
        assert!(w_zero.is_identity());
    }

    #[test]
    fn test_word_iterations_and_collections() {
        let terms = vec![Term::new(0, 1), Term::new(2, 3)];
        let w = Word::from_iter(terms.clone());
        assert_eq!(w.terms(), &terms[..]);

        let collected_refs: Vec<&Term> = w.iter().collect();
        assert_eq!(collected_refs, vec![&terms[0], &terms[1]]);

        let collected_owned: Vec<Term> = w.into_iter().collect();
        assert_eq!(collected_owned, terms);
    }
}
