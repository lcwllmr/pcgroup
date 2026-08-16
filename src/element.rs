//! Group element representation and arithmetic for polycyclic (PC) groups.
//!
//! An [`Element`] binds a canonical normal form [`Word`] to its associated [`Presentation`],
//! enabling intuitive group arithmetic (`*`, `pow`, `inverse`, `==`).

use std::fmt;
use std::ops::Mul;

use crate::word::Term;
use crate::{Collector, Presentation, Word};

/// A mathematical element of a polycyclic group.
///
/// Binds a normal-form [word](Word) to its [presentation](Presentation) so group operations
/// can be performed natively.
///
/// # Example: Arithmetic in S_3
/// ```
/// use pcgroup::{Builder, Element, Word};
///
/// // Presentation of S_3: < g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >
/// let pres = Builder::new(vec![2, 3])
///     .unwrap()
///     .add_commutator(1, 0, Word::from_term(1, 1))
///     .unwrap()
///     .build();
///
/// let g0 = Element::from_generator(0, 1, &pres);
/// let g1 = Element::from_generator(1, 1, &pres);
///
/// // In S_3: g1 * g0 = g0 * g1^2
/// let prod = &g1 * &g0;
/// assert_eq!(prod.to_string(), "g0 g1^2");
///
/// // g0 * g0 = 1 (order 2)
/// assert_eq!((&g0 * &g0).is_identity(), true);
///
/// // Inverses: (g0 * g1) * (g0 * g1)^-1 = 1
/// let elem = &g0 * &g1;
/// let elem_inv = elem.inverse();
/// assert_eq!((&elem * &elem_inv).is_identity(), true);
/// ```
#[derive(Debug, Clone)]
pub struct Element<'a> {
    pub word: Word,
    pub pres: &'a Presentation,
}

impl<'a> Element<'a> {
    /// Creates a new element. The word is assumed to already be in normal form.
    #[inline]
    pub fn new(word: Word, pres: &'a Presentation) -> Self {
        Self { word, pres }
    }

    /// Creates the identity element for the given presentation.
    #[inline]
    pub fn identity(pres: &'a Presentation) -> Self {
        Self {
            word: Word::identity(),
            pres,
        }
    }

    /// Creates an element consisting of a single generator power `gi^exponent`.
    #[inline]
    pub fn from_generator(gen_index: usize, exponent: u32, pres: &'a Presentation) -> Self {
        let mut collector = Collector::new(pres);
        collector.collect_generator(gen_index, exponent);
        Self {
            word: collector.into_word(),
            pres,
        }
    }

    /// Returns a reference to the underlying normal form [`Word`].
    #[inline]
    pub fn word(&self) -> &Word {
        &self.word
    }

    /// Returns a reference to the associated [`Presentation`].
    #[inline]
    pub fn presentation(&self) -> &'a Presentation {
        self.pres
    }

    /// Returns `true` if this element is the group identity.
    #[inline]
    pub fn is_identity(&self) -> bool {
        self.word.is_identity()
    }

    /// Returns `true` if the underlying word has no terms (equivalent to [`is_identity`](Self::is_identity)).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.word.is_empty()
    }

    /// Returns the number of terms in the underlying normal form word.
    #[inline]
    pub fn len(&self) -> usize {
        self.word.len()
    }

    /// Returns the leading term if non-identity.
    #[inline]
    pub fn leading_term(&self) -> Option<Term> {
        self.word.leading_term()
    }

    /// Returns the leading generator index if non-identity.
    #[inline]
    pub fn leading_generator(&self) -> Option<usize> {
        self.word.leading_generator()
    }

    /// Returns the leading generator exponent if non-identity.
    #[inline]
    pub fn leading_exponent(&self) -> Option<u32> {
        self.word.leading_exponent()
    }

    /// Computes `self^exp` using binary exponentiation.
    pub fn pow(&self, mut exp: u32) -> Self {
        if exp == 0 {
            return Self::identity(self.pres);
        }
        let mut base = self.clone();
        let mut result = Self::identity(self.pres);
        while exp > 0 {
            if exp % 2 == 1 {
                result = if result.is_identity() {
                    base.clone()
                } else {
                    &result * &base
                };
            }
            exp /= 2;
            if exp > 0 {
                base = &base * &base;
            }
        }
        result
    }

    /// Computes the inverse of the element using structural PC relations.
    ///
    /// For a word `x = g{i1}^{e1} ... g{ik}^{ek}`, the inverse is
    /// `x^-1 = (g{ik}^-1)^{ek} ... (g{i1}^-1)^{e1}`.
    ///
    /// Generator inverses are computed via `gi^-1 = gi^{pi - 1} * (gi^pi)^-1`,
    /// which terminates recursively because power tails only contain generators strictly greater than `i`.
    pub fn inverse(&self) -> Self {
        if self.is_identity() {
            return self.clone();
        }

        let mut result = Self::identity(self.pres);
        for term in self.word.iter() {
            let gen_inv = self.generator_inverse(term.gen_index);
            let term_inv = gen_inv.pow(term.exponent);
            result = term_inv * result;
        }
        result
    }

    /// Computes the inverse of a single standard generator `gi`.
    fn generator_inverse(&self, gen_idx: usize) -> Self {
        let p = self.pres.relative_order(gen_idx);
        let g_pow = if p > 1 {
            Self::new(Word::from_term(gen_idx, p - 1), self.pres)
        } else {
            Self::identity(self.pres)
        };

        let power_tail = self.pres.power(gen_idx);
        if power_tail.is_identity() {
            g_pow
        } else {
            let tail_elem = Self::new(power_tail.clone(), self.pres);
            let tail_inv = tail_elem.inverse();
            &g_pow * &tail_inv
        }
    }
}

impl<'a> fmt::Display for Element<'a> {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.word)
    }
}

impl<'a> PartialEq for Element<'a> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.word == other.word
    }
}

impl<'a> Eq for Element<'a> {}

impl<'a> Mul for Element<'a> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl<'a> Mul<&Element<'a>> for Element<'a> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &Element<'a>) -> Self::Output {
        &self * rhs
    }
}

impl<'a> Mul<Element<'a>> for &Element<'a> {
    type Output = Element<'a>;

    #[inline]
    fn mul(self, rhs: Element<'a>) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&Element<'a>> for &Element<'a> {
    type Output = Element<'a>;

    fn mul(self, rhs: &Element<'a>) -> Self::Output {
        let mut collector = Collector::new(self.pres);
        collector.collect(&self.word);
        collector.collect(&rhs.word);

        Element {
            word: collector.into_word(),
            pres: self.pres,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Builder;

    #[test]
    fn test_element_identity() {
        let pres = Builder::new(vec![2, 3]).unwrap().build();
        let id = Element::identity(&pres);

        assert!(id.is_identity());
        assert!(id.is_empty());
        assert_eq!(id.len(), 0);
        assert_eq!(id.to_string(), "1");
        assert_eq!(id.leading_term(), None);
        assert_eq!(id.leading_generator(), None);
        assert_eq!(id.leading_exponent(), None);
    }

    #[test]
    fn test_element_multiplication_and_powers() {
        let pres = Builder::new(vec![2, 3])
            .unwrap()
            .add_commutator(1, 0, Word::from_term(1, 1))
            .unwrap()
            .build();

        let g0 = Element::from_generator(0, 1, &pres);
        let g1 = Element::from_generator(1, 1, &pres);

        // In S_3: g1 * g0 = g0 * g1^2
        let prod = &g1 * &g0;
        assert_eq!(prod.to_string(), "g0 g1^2");

        // g0^2 = 1
        assert_eq!(g0.pow(2), Element::identity(&pres));
        assert_eq!(&g0 * &g0, Element::identity(&pres));

        // g1^3 = 1
        assert_eq!(g1.pow(3), Element::identity(&pres));

        // Associativity: (g1 * g0) * g1 == g1 * (g0 * g1)
        let left = (&g1 * &g0) * &g1;
        let right = &g1 * (&g0 * &g1);
        assert_eq!(left, right);
    }

    #[test]
    fn test_element_inverses_s3() {
        let pres = Builder::new(vec![2, 3])
            .unwrap()
            .add_commutator(1, 0, Word::from_term(1, 1))
            .unwrap()
            .build();

        let id = Element::identity(&pres);
        assert_eq!(id.inverse(), id);

        let g0 = Element::from_generator(0, 1, &pres);
        let g1 = Element::from_generator(1, 1, &pres);

        // g0^-1 = g0 (order 2)
        assert_eq!(g0.inverse(), g0);

        // g1^-1 = g1^2 (order 3)
        let g1_inv = g1.inverse();
        assert_eq!(g1_inv.to_string(), "g1^2");
        assert_eq!(&g1 * &g1_inv, id);
        assert_eq!(&g1_inv * &g1, id);

        // (g0 * g1)^-1
        let elem = &g0 * &g1;
        let elem_inv = elem.inverse();
        assert_eq!(&elem * &elem_inv, id);
        assert_eq!(&elem_inv * &elem, id);
    }

    #[test]
    fn test_element_inverses_c4() {
        // C_4: < g0, g1 | g0^2 = g1, g1^2 = 1, [g1, g0] = 1 >
        let pres = Builder::new(vec![2, 2])
            .unwrap()
            .add_power(0, Word::from_term(1, 1))
            .unwrap()
            .build();

        let id = Element::identity(&pres);
        let g0 = Element::from_generator(0, 1, &pres);

        // In C_4, generator g0 has order 4; g0^-1 = g0^3 = g0 * g1
        let g0_inv = g0.inverse();
        assert_eq!(g0_inv.to_string(), "g0 g1");
        assert_eq!(&g0 * &g0_inv, id);
        assert_eq!(&g0_inv * &g0, id);
    }
}
