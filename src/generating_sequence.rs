//! Polycyclic Generating Sequences (PGS / CGS) for subgroups of polycyclic groups.
//!
//! A [`GeneratingSequence`] maintains a sequence of polycyclic generator elements defining
//! a subgroup in echelon form with normalized leading exponents, providing fast membership
//! testing, subgroup closure (Sims' algorithm), order/index computation, and canonical reduction.

use std::fmt;

use crate::word::Term;
use crate::{Element, Presentation, Word};

/// A polycyclic generating sequence defining a subgroup of a polycyclic group.
///
/// Elements in a valid generating sequence are maintained in strict echelon form:
/// - Sorted in ascending order by their leading generator index (`i0 < i1 < ... < ik`).
/// - Each leading generator index is unique.
/// - Each element is normalized to have a leading exponent of `1`.
///
/// # Example: Subgroups of D_8
/// ```
/// use pcgroup::zoo::dihedral;
/// use pcgroup::{Element, GeneratingSequence};
///
/// // D_8 of order 8: reflection g0 (order 2) and rotation g1 (order 4)
/// let d8 = dihedral(4);
/// let s = Element::from_generator(0, 1, &d8);
/// let r = Element::from_generator(1, 1, &d8);
///
/// // Construct the cyclic rotation subgroup < r > = C_4
/// let rot_subgroup = GeneratingSequence::from_generators(&d8, &[r.clone()]);
/// assert_eq!(rot_subgroup.order(), 4);
/// assert_eq!(rot_subgroup.index(), 2);
///
/// // Membership testing
/// assert!(rot_subgroup.contains(&r));
/// assert!(!rot_subgroup.contains(&s));
/// ```
#[derive(Debug, Clone)]
pub struct GeneratingSequence<'a> {
    pres: &'a Presentation,
    elements: Vec<Element<'a>>,
}

impl<'a> GeneratingSequence<'a> {
    /// Initializes an empty generating sequence representing the trivial subgroup `{1}`.
    #[inline]
    pub fn trivial(pres: &'a Presentation) -> Self {
        Self {
            pres,
            elements: Vec::new(),
        }
    }

    /// Initializes an empty generating sequence (equivalent to [`trivial`](Self::trivial)).
    #[inline]
    pub fn empty(pres: &'a Presentation) -> Self {
        Self::trivial(pres)
    }

    /// Initializes a generating sequence representing the full group `G = < g0, ..., g{n-1} >`.
    pub fn full_group(pres: &'a Presentation) -> Self {
        let mut elements = Vec::with_capacity(pres.num_gens());
        for i in 0..pres.num_gens() {
            elements.push(Element::from_generator(i, 1, pres));
        }
        Self { pres, elements }
    }

    /// Constructs a closed polycyclic generating sequence from arbitrary generators using Sims' algorithm.
    ///
    /// The resulting sequence is closed under all power and commutator relations, and is
    /// tail-reduced into canonical generating sequence (CGS) form.
    pub fn from_generators(pres: &'a Presentation, gens: &[Element<'a>]) -> Self {
        let mut seq = Self::trivial(pres);
        for g in gens {
            seq.insert(g.clone());
        }

        let mut i = 0;
        while i < seq.elements.len() {
            // Power relation check: u_i^p_{lead(u_i)}
            let lead_gen = seq.elements[i].leading_generator().unwrap();
            let p = pres.relative_order(lead_gen);
            let pwr = seq.elements[i].pow(p);
            let rem_pwr = seq.sift(pwr);
            if !rem_pwr.is_identity() {
                seq.insert(rem_pwr);
                i = 0;
                continue;
            }

            // Commutator relations check: [u_j, u_i] = u_j^-1 * u_i^-1 * u_j * u_i for j > i
            let mut restart = false;
            for j in (i + 1)..seq.elements.len() {
                let u_i = &seq.elements[i];
                let u_j = &seq.elements[j];
                let comm = &u_j.inverse() * &(&u_i.inverse() * &(u_j * u_i));
                let rem_comm = seq.sift(comm);
                if !rem_comm.is_identity() {
                    seq.insert(rem_comm);
                    restart = true;
                    break;
                }
            }

            if restart {
                i = 0;
            } else {
                i += 1;
            }
        }

        seq.canonicalize();
        seq
    }

    /// Returns a reference to the associated [`Presentation`].
    #[inline]
    pub fn presentation(&self) -> &'a Presentation {
        self.pres
    }

    /// Returns a slice over the basis [`Element`]s in this sequence.
    #[inline]
    pub fn elements(&self) -> &[Element<'a>] {
        &self.elements
    }

    /// Returns `true` if this generating sequence represents the trivial subgroup `{1}`.
    #[inline]
    pub fn is_trivial(&self) -> bool {
        self.elements.is_empty()
    }

    /// Returns `true` if the sequence has no elements (equivalent to [`is_trivial`](Self::is_trivial)).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Returns the number of polycyclic generators in this sequence.
    #[inline]
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// Returns an iterator over references to the basis elements.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, Element<'a>> {
        self.elements.iter()
    }

    /// Computes the order of this subgroup `|H| = prod p_{lead(u_i)}`.
    #[inline]
    pub fn order(&self) -> u128 {
        self.elements
            .iter()
            .map(|e| self.pres.relative_order(e.leading_generator().unwrap()) as u128)
            .product()
    }

    /// Computes the index `[G : H] = |G| / |H|` of this subgroup in the parent group.
    #[inline]
    pub fn index(&self) -> u128 {
        self.pres.order() / self.order()
    }

    /// Returns `true` if element `x` is contained in this subgroup.
    #[inline]
    pub fn contains(&self, x: &Element<'a>) -> bool {
        self.sift(x.clone()).is_identity()
    }

    /// Sifts an element `x` through the generating sequence, eliminating matching leading terms.
    ///
    /// Returns the remainder (which is `1` if `x` belongs to the subgroup).
    pub fn sift(&self, mut x: Element<'a>) -> Element<'a> {
        while let Some(term) = x.leading_term() {
            if let Some(h) = self
                .elements
                .iter()
                .find(|e| e.leading_generator() == Some(term.gen_index))
            {
                let h_inv = h.inverse().pow(term.exponent);
                x = h_inv * x;
            } else {
                break;
            }
        }
        x
    }

    /// Inserts an element `x` into the sequence while maintaining normalized echelon form.
    ///
    /// Returns `true` if the subgroup grew, or `false` if `x` was already in the span.
    pub fn insert(&mut self, x: Element<'a>) -> bool {
        let rem = self.sift(x);
        if rem.is_identity() {
            return false;
        }

        let lead_term = rem.leading_term().unwrap();
        let p = self.pres.relative_order(lead_term.gen_index);
        let inv = crate::util::mod_inverse(lead_term.exponent, p)
            .expect("leading exponent must be coprime to prime relative order");
        let norm_rem = rem.pow(inv);

        let idx = self
            .elements
            .partition_point(|e| e.leading_generator().unwrap() < lead_term.gen_index);
        self.elements.insert(idx, norm_rem);
        true
    }

    /// Reduces trailing terms of all basis elements to obtain a unique Canonical Generating Sequence (CGS).
    pub fn canonicalize(&mut self) {
        let n = self.elements.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let target_gen = self.elements[j].leading_generator().unwrap();
                if let Some(exp) = self.elements[i]
                    .word()
                    .iter()
                    .find(|t| t.gen_index == target_gen)
                    .map(|t| t.exponent)
                {
                    let u_j_inv = self.elements[j].inverse().pow(exp);
                    self.elements[i] = &self.elements[i] * &u_j_inv;
                }
            }
        }
    }

    /// Expresses an element `x` as a [`Word`] over the basis indices `0..self.len()`.
    ///
    /// Returns `Some(word)` if `x` is in this subgroup, or `None` if `x` is not in this subgroup.
    pub fn express_in_basis(&self, mut x: Element<'a>) -> Option<Word> {
        let mut terms = Vec::with_capacity(self.elements.len());
        for (pos, b) in self.elements.iter().enumerate() {
            let b_lead = b.leading_generator().unwrap();
            if let Some(x_lead) = x.leading_generator() {
                if x_lead == b_lead {
                    let exp = x.leading_exponent().unwrap();
                    terms.push(Term::new(pos, exp));
                    let b_inv = b.inverse().pow(exp);
                    x = b_inv * x;
                } else if x_lead < b_lead {
                    return None;
                }
            }
        }
        if x.is_identity() {
            Some(Word::new(terms))
        } else {
            None
        }
    }
}

impl<'a> fmt::Display for GeneratingSequence<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.elements.is_empty() {
            write!(f, "< 1 >")
        } else {
            write!(f, "< ")?;
            for (i, elem) in self.elements.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", elem)?;
            }
            write!(f, " >")
        }
    }
}

impl<'a> PartialEq for GeneratingSequence<'a> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.pres == other.pres && self.elements == other.elements
    }
}

impl<'a> Eq for GeneratingSequence<'a> {}

impl<'a> IntoIterator for GeneratingSequence<'a> {
    type Item = Element<'a>;
    type IntoIter = std::vec::IntoIter<Element<'a>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.elements.into_iter()
    }
}

impl<'a> IntoIterator for &'a GeneratingSequence<'a> {
    type Item = &'a Element<'a>;
    type IntoIter = std::slice::Iter<'a, Element<'a>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.elements.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zoo::{cyclic, dihedral, quaternion};

    #[test]
    fn test_trivial_and_full_group() {
        let pres = cyclic(6);
        let trivial = GeneratingSequence::trivial(&pres);
        assert!(trivial.is_trivial());
        assert_eq!(trivial.len(), 0);
        assert_eq!(trivial.order(), 1);
        assert_eq!(trivial.index(), 6);
        assert_eq!(trivial.to_string(), "< 1 >");

        let full = GeneratingSequence::full_group(&pres);
        assert!(!full.is_trivial());
        assert_eq!(full.len(), 2);
        assert_eq!(full.order(), 6);
        assert_eq!(full.index(), 1);
    }

    #[test]
    fn test_sift_and_normalization() {
        let pres = cyclic(5);
        let mut seq = GeneratingSequence::trivial(&pres);

        // Insert g0^2 in C_5: should normalize to g0^1 since 2^-1 = 3 mod 5 and (g0^2)^3 = g0^6 = g0^1
        let elem2 = Element::from_generator(0, 2, &pres);
        assert!(seq.insert(elem2));
        assert_eq!(seq.len(), 1);
        assert_eq!(seq.elements()[0].to_string(), "g0");
        assert_eq!(seq.order(), 5);

        // Sifting g0^4 should yield identity
        let elem4 = Element::from_generator(0, 4, &pres);
        assert!(seq.sift(elem4).is_identity());
        assert!(seq.contains(&Element::from_generator(0, 3, &pres)));
    }

    #[test]
    fn test_subgroup_closure_dihedral() {
        // D_8 of order 8: g0 is reflection s (order 2), g1 is rotation r (order 4)
        let pres = dihedral(4);
        let s = Element::from_generator(0, 1, &pres); // reflection of order 2
        let r = Element::from_generator(1, 1, &pres); // rotation of order 4

        // Subgroup generated by reflection s: order 2
        let sub_s = GeneratingSequence::from_generators(&pres, &[s.clone()]);
        assert_eq!(sub_s.order(), 2);
        assert_eq!(sub_s.index(), 4);
        assert!(sub_s.contains(&s));
        assert!(!sub_s.contains(&r));

        // Subgroup generated by rotation r: order 4
        let sub_r = GeneratingSequence::from_generators(&pres, &[r.clone()]);
        assert_eq!(sub_r.order(), 4);
        assert_eq!(sub_r.index(), 2);
        assert!(sub_r.contains(&r));
        assert!(!sub_r.contains(&s));

        // Subgroup generated by r and s: full group of order 8
        let sub_full = GeneratingSequence::from_generators(&pres, &[r.clone(), s.clone()]);
        assert_eq!(sub_full.order(), 8);
        assert_eq!(sub_full.index(), 1);
    }

    #[test]
    fn test_subgroup_closure_quaternion() {
        // Q_8 of order 8: generators y = g0, x = g1
        let pres = quaternion(2);
        let y_elem = Element::from_generator(0, 1, &pres);
        let x_elem = Element::from_generator(1, 1, &pres);

        // Cyclic subgroup <x> of order 4
        let sub_x = GeneratingSequence::from_generators(&pres, &[x_elem.clone()]);
        assert_eq!(sub_x.order(), 4);
        assert_eq!(sub_x.index(), 2);
        assert!(sub_x.contains(&x_elem));
        assert!(sub_x.contains(&x_elem.pow(2)));
        assert!(!sub_x.contains(&y_elem));

        // Subgroup generated by <x> and <y>: full Q_8 of order 8
        let sub_q8 = GeneratingSequence::from_generators(&pres, &[x_elem, y_elem]);
        assert_eq!(sub_q8.order(), 8);
        assert_eq!(sub_q8.index(), 1);
    }

    #[test]
    fn test_canonicalize_and_equality() {
        let pres = cyclic(6);
        let g0 = Element::from_generator(0, 1, &pres);
        let g1 = Element::from_generator(1, 1, &pres);

        // Generate full group using two different generator sets
        let seq1 = GeneratingSequence::from_generators(&pres, &[&g0 * &g1, g1.clone()]);
        let seq2 = GeneratingSequence::full_group(&pres);
        assert_eq!(seq1, seq2);
    }

    #[test]
    fn test_express_in_basis() {
        let pres = dihedral(4);
        let s = Element::from_generator(0, 1, &pres);
        let r = Element::from_generator(1, 1, &pres);
        let full = GeneratingSequence::full_group(&pres);

        let elem = &s * &r;
        let word = full.express_in_basis(elem.clone()).unwrap();
        assert_eq!(word, elem.word);

        // Subgroup <s>
        let sub_s = GeneratingSequence::from_generators(&pres, &[s.clone()]);
        assert_eq!(sub_s.express_in_basis(s).unwrap().to_string(), "g0");
        assert!(sub_s.express_in_basis(r).is_none());
    }
}
