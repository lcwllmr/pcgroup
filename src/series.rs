//! Subgroup series and structural property checks for polycyclic groups.
//!
//! This module provides algorithms for computing descending subgroup series such as the
//! [lower central series](lower_central_series), computing [commutator subgroups](commutator_subgroup),
//! and testing structural properties including [abelianity](is_abelian) and [nilpotency](nilpotency_class).

use crate::{Element, GeneratingSequence};

/// Computes the commutator subgroup `[h_group, k_group]` closed under `group_gens`.
///
/// # Algorithm
/// 1. Initialize an empty [`GeneratingSequence`] `c`.
/// 2. For every generator `h` in `h_group` and `k` in `k_group`, compute the commutator
///    `[h, k] = h^-1 * k^-1 * h * k`.
/// 3. Sift and insert each commutator into `c`.
/// 4. Return the normal closure `c.normal_closure_in(group_gens)`.
///
/// # Examples
/// ```
/// use pcgroup::zoo::{dihedral, cyclic};
/// use pcgroup::{Element, GeneratingSequence, commutator_subgroup};
///
/// // D_8 of order 8
/// let d8 = dihedral(4);
/// let full = GeneratingSequence::full_group(&d8);
/// let full_gens: Vec<_> = (0..d8.num_gens())
///     .map(|i| Element::from_generator(i, 1, &d8))
///     .collect();
///
/// // [D_8, D_8] is the center Z(D_8) of order 2
/// let derived = commutator_subgroup(&full, &full, &full_gens);
/// assert_eq!(derived.order(), 2);
/// ```
pub fn commutator_subgroup<'a>(
    h_group: &GeneratingSequence<'a>,
    k_group: &GeneratingSequence<'a>,
    group_gens: &[Element<'a>],
) -> GeneratingSequence<'a> {
    let mut c = GeneratingSequence::empty(h_group.presentation());

    for h in h_group.elements() {
        let h_inv = h.inverse();
        for k in k_group.elements() {
            let k_inv = k.inverse();
            // [h, k] = h^-1 * k^-1 * h * k
            let comm = &h_inv * &k_inv * h * k;
            c.insert(comm);
        }
    }

    c.normal_closure_in(group_gens)
}

/// Checks if a group or subgroup is abelian.
///
/// A group `H` is abelian if and only if its derived (commutator) subgroup `[H, H]` is trivial (`{1}`).
///
/// # Examples
/// ```
/// use pcgroup::zoo::{cyclic, dihedral};
/// use pcgroup::{Element, GeneratingSequence, is_abelian};
///
/// let c6 = cyclic(6);
/// let full_c6 = GeneratingSequence::full_group(&c6);
/// let c6_gens: Vec<_> = (0..c6.num_gens()).map(|i| Element::from_generator(i, 1, &c6)).collect();
/// assert!(is_abelian(&full_c6, &c6_gens));
///
/// let d6 = dihedral(3); // S_3 of order 6 (non-abelian)
/// let full_d6 = GeneratingSequence::full_group(&d6);
/// let d6_gens: Vec<_> = (0..d6.num_gens()).map(|i| Element::from_generator(i, 1, &d6)).collect();
/// assert!(!is_abelian(&full_d6, &d6_gens));
/// ```
#[inline]
pub fn is_abelian<'a>(group: &GeneratingSequence<'a>, group_gens: &[Element<'a>]) -> bool {
    commutator_subgroup(group, group, group_gens).is_trivial()
}

/// Computes a descending series of subgroups starting from `start` and iteratively applying `step`.
///
/// Iteration stops when the current subgroup reaches the trivial group `{1}` or stabilizes (`step(curr) == curr`).
pub(crate) fn compute_descending_series<'a, F>(
    start: GeneratingSequence<'a>,
    mut step: F,
) -> Vec<GeneratingSequence<'a>>
where
    F: FnMut(&GeneratingSequence<'a>) -> GeneratingSequence<'a>,
{
    let mut series = vec![start];
    loop {
        let curr = series.last().unwrap();
        if curr.is_trivial() {
            break;
        }

        let next = step(curr);
        if next == *curr {
            break;
        }

        let is_triv = next.is_trivial();
        series.push(next);
        if is_triv {
            break;
        }
    }
    series
}

/// Computes the lower central series of a polycyclic group or subgroup.
///
/// The lower central series is defined as:
/// - `gamma_1 = group`
/// - `gamma_{i+1} = [gamma_i, group]` (closed under `group_gens`)
///
/// Iteration stops when `gamma_{i+1}` stabilizes or becomes trivial.
///
/// # Examples
/// ```
/// use pcgroup::zoo::dihedral;
/// use pcgroup::{Element, GeneratingSequence, lower_central_series};
///
/// // D_8 of order 8: gamma_1 = D_8 (order 8), gamma_2 = Z(D_8) (order 2), gamma_3 = {1} (order 1)
/// let d8 = dihedral(4);
/// let full = GeneratingSequence::full_group(&d8);
/// let full_gens: Vec<_> = (0..d8.num_gens()).map(|i| Element::from_generator(i, 1, &d8)).collect();
///
/// let lcs = lower_central_series(&full, &full_gens);
/// assert_eq!(lcs.len(), 3);
/// assert_eq!(lcs[0].order(), 8);
/// assert_eq!(lcs[1].order(), 2);
/// assert_eq!(lcs[2].order(), 1);
/// ```
pub fn lower_central_series<'a>(
    group: &GeneratingSequence<'a>,
    group_gens: &[Element<'a>],
) -> Vec<GeneratingSequence<'a>> {
    compute_descending_series(group.clone(), |curr| {
        commutator_subgroup(curr, group, group_gens)
    })
}

/// Computes the nilpotency class of a polycyclic group or subgroup.
///
/// Returns `Some(c)` where `c` is the nilpotency class if the group is nilpotent:
/// - `c = 0` if `group` is the trivial group `{1}`.
/// - `c = 1` if `group` is non-trivial and abelian (`[G, G] = {1}`).
/// - `c >= 2` if the lower central series terminates at `{1}` in `c + 1` terms (`gamma_{c+1} = {1}` with `gamma_c != {1}`).
///
/// Returns `None` if the group is not nilpotent (i.e. the lower central series stabilizes at a non-trivial subgroup).
///
/// # Examples
/// ```
/// use pcgroup::zoo::{abelian, dihedral, quaternion};
/// use pcgroup::{Element, GeneratingSequence, nilpotency_class};
///
/// // Klein 4-group: abelian, nilpotency class 1
/// let v4 = abelian(&[2, 2]);
/// let full_v4 = GeneratingSequence::full_group(&v4);
/// let v4_gens: Vec<_> = (0..v4.num_gens()).map(|i| Element::from_generator(i, 1, &v4)).collect();
/// assert_eq!(nilpotency_class(&full_v4, &v4_gens), Some(1));
///
/// // Q_8: class 2 nilpotent
/// let q8 = quaternion(2);
/// let full_q8 = GeneratingSequence::full_group(&q8);
/// let q8_gens: Vec<_> = (0..q8.num_gens()).map(|i| Element::from_generator(i, 1, &q8)).collect();
/// assert_eq!(nilpotency_class(&full_q8, &q8_gens), Some(2));
///
/// // S_3 = D_6: non-nilpotent (stabilizes at A_3)
/// let s3 = dihedral(3);
/// let full_s3 = GeneratingSequence::full_group(&s3);
/// let s3_gens: Vec<_> = (0..s3.num_gens()).map(|i| Element::from_generator(i, 1, &s3)).collect();
/// assert_eq!(nilpotency_class(&full_s3, &s3_gens), None);
/// ```
pub fn nilpotency_class<'a>(
    group: &GeneratingSequence<'a>,
    group_gens: &[Element<'a>],
) -> Option<usize> {
    let lcs = lower_central_series(group, group_gens);
    if lcs.last().is_some_and(|g| g.is_trivial()) {
        Some(lcs.len() - 1)
    } else {
        None
    }
}

/// Checks whether a polycyclic group or subgroup is nilpotent.
///
/// Returns `true` if the nilpotency class is finite, or `false` otherwise.
///
/// # Examples
/// ```
/// use pcgroup::zoo::{dihedral, quaternion};
/// use pcgroup::{Element, GeneratingSequence, is_nilpotent};
///
/// let d8 = dihedral(4);
/// let full_d8 = GeneratingSequence::full_group(&d8);
/// let d8_gens: Vec<_> = (0..d8.num_gens()).map(|i| Element::from_generator(i, 1, &d8)).collect();
/// assert!(is_nilpotent(&full_d8, &d8_gens));
///
/// let s3 = dihedral(3);
/// let full_s3 = GeneratingSequence::full_group(&s3);
/// let s3_gens: Vec<_> = (0..s3.num_gens()).map(|i| Element::from_generator(i, 1, &s3)).collect();
/// assert!(!is_nilpotent(&full_s3, &s3_gens));
/// ```
#[inline]
pub fn is_nilpotent<'a>(group: &GeneratingSequence<'a>, group_gens: &[Element<'a>]) -> bool {
    nilpotency_class(group, group_gens).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zoo::{abelian, cyclic, dihedral, quaternion};

    #[test]
    fn test_trivial_group_properties() {
        let pres = cyclic(1);
        let triv = GeneratingSequence::trivial(&pres);
        let gens = vec![];

        assert!(is_abelian(&triv, &gens));
        assert!(is_nilpotent(&triv, &gens));
        assert_eq!(nilpotency_class(&triv, &gens), Some(0));

        let lcs = lower_central_series(&triv, &gens);
        assert_eq!(lcs.len(), 1);
        assert!(lcs[0].is_trivial());
    }

    #[test]
    fn test_cyclic_and_abelian_series() {
        let pres = cyclic(12);
        let full = GeneratingSequence::full_group(&pres);
        let gens: Vec<_> = (0..pres.num_gens())
            .map(|i| Element::from_generator(i, 1, &pres))
            .collect();

        assert!(is_abelian(&full, &gens));
        assert_eq!(nilpotency_class(&full, &gens), Some(1));

        let lcs = lower_central_series(&full, &gens);
        assert_eq!(lcs.len(), 2);
        assert_eq!(lcs[0].order(), 12);
        assert_eq!(lcs[1].order(), 1);

        let v4 = abelian(&[2, 2]);
        let full_v4 = GeneratingSequence::full_group(&v4);
        let v4_gens: Vec<_> = (0..v4.num_gens())
            .map(|i| Element::from_generator(i, 1, &v4))
            .collect();

        assert!(is_abelian(&full_v4, &v4_gens));
        assert_eq!(nilpotency_class(&full_v4, &v4_gens), Some(1));
    }

    #[test]
    fn test_dihedral_series_and_nilpotency() {
        // D_8 of order 8 (nilpotent of class 2)
        let d8 = dihedral(4);
        let full_d8 = GeneratingSequence::full_group(&d8);
        let d8_gens: Vec<_> = (0..d8.num_gens())
            .map(|i| Element::from_generator(i, 1, &d8))
            .collect();

        assert!(!is_abelian(&full_d8, &d8_gens));
        assert!(is_nilpotent(&full_d8, &d8_gens));
        assert_eq!(nilpotency_class(&full_d8, &d8_gens), Some(2));

        let lcs_d8 = lower_central_series(&full_d8, &d8_gens);
        assert_eq!(lcs_d8.len(), 3);
        assert_eq!(lcs_d8[0].order(), 8);
        assert_eq!(lcs_d8[1].order(), 2);
        assert_eq!(lcs_d8[2].order(), 1);

        // S_3 = D_6 of order 6 (non-nilpotent)
        let s3 = dihedral(3);
        let full_s3 = GeneratingSequence::full_group(&s3);
        let s3_gens: Vec<_> = (0..s3.num_gens())
            .map(|i| Element::from_generator(i, 1, &s3))
            .collect();

        assert!(!is_abelian(&full_s3, &s3_gens));
        assert!(!is_nilpotent(&full_s3, &s3_gens));
        assert_eq!(nilpotency_class(&full_s3, &s3_gens), None);

        let lcs_s3 = lower_central_series(&full_s3, &s3_gens);
        // gamma_1 = S_3 (order 6), gamma_2 = A_3 (order 3), gamma_3 = [A_3, S_3] = A_3 (stabilizes)
        assert_eq!(lcs_s3.len(), 2);
        assert_eq!(lcs_s3[0].order(), 6);
        assert_eq!(lcs_s3[1].order(), 3);
    }

    #[test]
    fn test_quaternion_series_and_nilpotency() {
        // Q_8 of order 8 (nilpotent of class 2)
        let q8 = quaternion(2);
        let full_q8 = GeneratingSequence::full_group(&q8);
        let q8_gens: Vec<_> = (0..q8.num_gens())
            .map(|i| Element::from_generator(i, 1, &q8))
            .collect();

        assert!(!is_abelian(&full_q8, &q8_gens));
        assert!(is_nilpotent(&full_q8, &q8_gens));
        assert_eq!(nilpotency_class(&full_q8, &q8_gens), Some(2));

        let lcs_q8 = lower_central_series(&full_q8, &q8_gens);
        assert_eq!(lcs_q8.len(), 3);
        assert_eq!(lcs_q8[0].order(), 8);
        assert_eq!(lcs_q8[1].order(), 2);
        assert_eq!(lcs_q8[2].order(), 1);

        // Q_12 of order 12 (dicyclic with n=3, non-nilpotent)
        let q12 = quaternion(3);
        let full_q12 = GeneratingSequence::full_group(&q12);
        let q12_gens: Vec<_> = (0..q12.num_gens())
            .map(|i| Element::from_generator(i, 1, &q12))
            .collect();

        assert!(!is_abelian(&full_q12, &q12_gens));
        assert!(!is_nilpotent(&full_q12, &q12_gens));
        assert_eq!(nilpotency_class(&full_q12, &q12_gens), None);
    }
}
