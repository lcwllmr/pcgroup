//! Extensions of polycyclic groups: split extensions (semidirect products) and central extensions.
//!
//! A group extension of a normal subgroup `N` by a quotient group `Q` is a group `G` fitting into
//! a short exact sequence:
//!
//! `1 → N → G → Q → 1`
//!
//! Under power-commutator (PC) presentations, generators of `G` are the concatenated generators
//! `g0, ..., g{m-1}` (from `Q`) followed by `g_m, ..., g_{m+k-1}` (from `N`), where `m = |Q.gens|`
//! and `k = |N.gens|`.

use std::fmt;

use crate::series::is_abelian;
use crate::word::Term;
use crate::{Builder, ConsistencyError, Element, Presentation, Word, verify_consistency};

/// Errors that can occur during the construction of a group extension.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExtensionError {
    /// Kernel group `N` is not abelian in a central extension.
    NonAbelianKernel,

    /// Action or tails slice dimension does not match the expected generator or relation count.
    DimensionMismatch { expected: usize, actual: usize },

    /// A relation word, tail word, or action word contains a generator index out of bounds for `N`.
    GeneratorOutOfBounds { index: usize, max: usize },

    /// The resulting PC presentation is inconsistent (fails confluence / critical overlap checks).
    Inconsistent(ConsistencyError),
}

impl fmt::Display for ExtensionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonAbelianKernel => {
                write!(f, "Central extension requires an abelian kernel group N")
            }
            Self::DimensionMismatch { expected, actual } => {
                write!(
                    f,
                    "Dimension mismatch: expected {expected} elements, got {actual}"
                )
            }
            Self::GeneratorOutOfBounds { index, max } => {
                write!(
                    f,
                    "Generator index {index} is out of bounds (maximum valid index is {max})"
                )
            }
            Self::Inconsistent(err) => {
                write!(f, "Extension presentation is inconsistent: {err}")
            }
        }
    }
}

impl std::error::Error for ExtensionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Inconsistent(err) => Some(err),
            _ => None,
        }
    }
}

impl From<ConsistencyError> for ExtensionError {
    fn from(err: ConsistencyError) -> Self {
        Self::Inconsistent(err)
    }
}

/// Shifts all generator indices in a normal form word by `shift`.
#[inline]
fn shift_word(word: &Word, shift: usize) -> Word {
    if shift == 0 || word.is_identity() {
        return word.clone();
    }
    let terms = word
        .iter()
        .map(|t| Term::new(t.gen_index + shift, t.exponent))
        .collect();
    Word::new(terms)
}

/// Evaluates an endomorphism of presentation `n` on a `word` in `n`, given generator images `images`.
fn eval_hom(word: &Word, images: &[Word], n: &Presentation) -> Word {
    let mut res = Element::identity(n);
    for term in word.iter() {
        let img = Element::new(images[term.gen_index].clone(), n);
        res = res * img.pow(term.exponent);
    }
    res.word
}

/// Computes the inverse image of generator `n_i` under an automorphism of `N`.
fn invert_automorphism_on_generator(
    gen_idx: usize,
    images: &[Word],
    n: &Presentation,
) -> Result<Word, ExtensionError> {
    let initial = Word::from_term(gen_idx, 1);
    let mut current = initial.clone();
    let max_iter = n.order() as usize;

    for _ in 0..max_iter.max(1) {
        let next = eval_hom(&current, images, n);
        if next == initial {
            return Ok(current);
        }
        current = next;
    }
    Err(ExtensionError::Inconsistent(
        ConsistencyError::PowerAndItself { i: gen_idx },
    ))
}

/// Constructs a split extension (semidirect product) `N ⋊_ϕ Q` of a normal subgroup `N` by a quotient group `Q`.
///
/// The action `action` specifies the automorphism `ϕ(q_j) ∈ Aut(N)` for each generator `q_j` of `Q` (`0 <= j < m`),
/// given as a list of images `ϕ(q_j)(n_i) = q_j * n_i * q_j^-1 ∈ N` for each generator `n_i` of `N` (`0 <= i < k`).
///
/// If `action` is empty, the trivial action is used, yielding the direct product `N × Q`.
///
/// # Examples
/// ```
/// use pcgroup::split_extension;
/// use pcgroup::zoo::cyclic;
/// use pcgroup::{Word, is_abelian, verify_consistency};
///
/// // Construct S_3 ≅ C_3 ⋊ C_2
/// let c3 = cyclic(3);
/// let c2 = cyclic(2);
///
/// // Inversion action: q0 * n0 * q0^-1 = n0^2
/// let action = vec![vec![Word::from_term(0, 2)]];
/// let s3 = split_extension(&c3, &c2, &action).unwrap();
///
/// assert_eq!(s3.order(), 6);
/// assert_eq!(s3.num_gens(), 2);
/// assert_eq!(verify_consistency(&s3), Ok(()));
/// assert!(!is_abelian(&s3));
/// ```
///
/// # Errors
/// Returns [`ExtensionError::DimensionMismatch`] if `action` length does not match `Q.num_gens()` or if
/// any inner list length does not match `N.num_gens()`.
/// Returns [`ExtensionError::GeneratorOutOfBounds`] if any action word references generator indices `>= N.num_gens()`.
/// Returns [`ExtensionError::Inconsistent`] if the resulting PC presentation is not confluent.
pub fn split_extension(
    n: &Presentation,
    q: &Presentation,
    action: &[Vec<Word>],
) -> Result<Presentation, ExtensionError> {
    let m = q.num_gens();
    let k = n.num_gens();

    if !action.is_empty() {
        if action.len() != m {
            return Err(ExtensionError::DimensionMismatch {
                expected: m,
                actual: action.len(),
            });
        }
        for q_act in action {
            if q_act.len() != k {
                return Err(ExtensionError::DimensionMismatch {
                    expected: k,
                    actual: q_act.len(),
                });
            }
            for word in q_act {
                for term in word.iter() {
                    if term.gen_index >= k {
                        return Err(ExtensionError::GeneratorOutOfBounds {
                            index: term.gen_index,
                            max: k.saturating_sub(1),
                        });
                    }
                }
            }
        }
    }

    let mut rel_orders = q.relative_orders().to_vec();
    rel_orders.extend_from_slice(n.relative_orders());

    let mut builder = Builder::new(rel_orders).expect("relative orders are primes");

    // 1. Power relations of Q (unchanged)
    for j in 0..m {
        builder = builder
            .add_power(j, q.power(j).clone())
            .expect("valid Q power");
    }

    // 2. Commutator relations of Q (unchanged)
    for j1 in 1..m {
        for j2 in 0..j1 {
            builder = builder
                .add_commutator(j1, j2, q.commutator(j1, j2).clone())
                .expect("valid Q commutator");
        }
    }

    // 3. Power relations of N (shifted by +m)
    for i in 0..k {
        let shifted = shift_word(n.power(i), m);
        builder = builder
            .add_power(m + i, shifted)
            .expect("valid N power relation");
    }

    // 4. Commutator relations of N (shifted by +m)
    for i1 in 1..k {
        for i2 in 0..i1 {
            let shifted = shift_word(n.commutator(i1, i2), m);
            builder = builder
                .add_commutator(m + i1, m + i2, shifted)
                .expect("valid N commutator relation");
        }
    }

    // 5. Cross action relations: [g_{m+i}, g_j] = n_i^-1 * ϕ(q_j^-1)(n_i)
    if !action.is_empty() {
        for (j, q_act) in action.iter().enumerate() {
            for i in 0..k {
                let inv_img = invert_automorphism_on_generator(i, q_act, n)?;
                let n_i_elem = Element::from_generator(i, 1, n);
                let inv_img_elem = Element::new(inv_img, n);
                let comm_elem = n_i_elem.inverse() * inv_img_elem;
                let shifted_comm = shift_word(&comm_elem.word, m);
                builder = builder
                    .add_commutator(m + i, j, shifted_comm)
                    .expect("valid cross commutator");
            }
        }
    }

    let pres = builder.build();
    verify_consistency(&pres)?;
    Ok(pres)
}

/// Constructs a central extension of an abelian normal subgroup `N` by a quotient group `Q` defined by relation tails.
///
/// The central extension fits into `1 → N → G → Q → 1` with `N ⊆ Z(G)`.
/// Power and commutator relations of `Q` hold modulo `N`, modified by appending elements of `N` ("tails"):
/// - `g_j^{p_j} = w_j * t_pow(j)` for `0 <= j < m`
/// - `[g_{j1}, g_{j2}] = w_{j1, j2} * t_comm(j1, j2)` for `0 <= j2 < j1 < m`
///
/// # Arguments
/// - `n`: The abelian normal subgroup `N`.
/// - `q`: The quotient group `Q`.
/// - `power_tails`: Tails in `N` for each power relation `q_j^{p_j}` of `Q` (`0 <= j < m`).
/// - `commutator_tails`: Tails in `N` for each commutator relation `[q_{j1}, q_{j2}]` of `Q` (`0 <= j2 < j1 < m`),
///   indexed according to standard lower-triangular layout [`Presentation::commutator_index`].
///
/// If `power_tails` or `commutator_tails` is empty, trivial tails (the identity word) are used.
///
/// # Examples
/// ```
/// use pcgroup::central_extension;
/// use pcgroup::zoo::{abelian, cyclic};
/// use pcgroup::{Word, is_abelian, is_nilpotent, nilpotency_class, verify_consistency};
///
/// // Construct Q_8 as a central extension of C_2 by V_4 = C_2 x C_2
/// let c2 = cyclic(2);
/// let v4 = abelian(&[2, 2]);
///
/// // Non-trivial power and commutator tails in C_2
/// let pow_tails = vec![Word::from_term(0, 1), Word::from_term(0, 1)];
/// let comm_tails = vec![Word::from_term(0, 1)];
///
/// let q8 = central_extension(&c2, &v4, &pow_tails, &comm_tails).unwrap();
/// assert_eq!(q8.order(), 8);
/// assert_eq!(verify_consistency(&q8), Ok(()));
/// assert!(!is_abelian(&q8));
/// assert!(is_nilpotent(&q8));
/// assert_eq!(nilpotency_class(&q8), Some(2));
/// ```
///
/// # Errors
/// Returns [`ExtensionError::NonAbelianKernel`] if `N` is not abelian.
/// Returns [`ExtensionError::DimensionMismatch`] if `power_tails` length does not match `Q.num_gens()` or
/// `commutator_tails` length does not match `Presentation::num_commutators(Q.num_gens())`.
/// Returns [`ExtensionError::GeneratorOutOfBounds`] if any tail word references generator indices `>= N.num_gens()`.
/// Returns [`ExtensionError::Inconsistent`] if the resulting PC presentation is not confluent.
pub fn central_extension(
    n: &Presentation,
    q: &Presentation,
    power_tails: &[Word],
    commutator_tails: &[Word],
) -> Result<Presentation, ExtensionError> {
    if !is_abelian(n) {
        return Err(ExtensionError::NonAbelianKernel);
    }

    let m = q.num_gens();
    let k = n.num_gens();
    let num_q_comms = Presentation::num_commutators(m);

    if !power_tails.is_empty() {
        if power_tails.len() != m {
            return Err(ExtensionError::DimensionMismatch {
                expected: m,
                actual: power_tails.len(),
            });
        }
        for word in power_tails {
            for term in word.iter() {
                if term.gen_index >= k {
                    return Err(ExtensionError::GeneratorOutOfBounds {
                        index: term.gen_index,
                        max: k.saturating_sub(1),
                    });
                }
            }
        }
    }

    if !commutator_tails.is_empty() {
        if commutator_tails.len() != num_q_comms {
            return Err(ExtensionError::DimensionMismatch {
                expected: num_q_comms,
                actual: commutator_tails.len(),
            });
        }
        for word in commutator_tails {
            for term in word.iter() {
                if term.gen_index >= k {
                    return Err(ExtensionError::GeneratorOutOfBounds {
                        index: term.gen_index,
                        max: k.saturating_sub(1),
                    });
                }
            }
        }
    }

    let mut rel_orders = q.relative_orders().to_vec();
    rel_orders.extend_from_slice(n.relative_orders());

    let mut builder = Builder::new(rel_orders).expect("relative orders are primes");

    // 1. Power relations of Q with appended tails in N
    for j in 0..m {
        let mut terms = q.power(j).terms.clone();
        if let Some(tail) = power_tails.get(j) {
            terms.extend(shift_word(tail, m).terms);
        }
        builder = builder
            .add_power(j, Word::new(terms))
            .expect("valid modified Q power");
    }

    // 2. Commutator relations of Q with appended tails in N
    for j1 in 1..m {
        for j2 in 0..j1 {
            let idx = Presentation::commutator_index(j1, j2);
            let mut terms = q.commutator(j1, j2).terms.clone();
            if let Some(tail) = commutator_tails.get(idx) {
                terms.extend(shift_word(tail, m).terms);
            }
            builder = builder
                .add_commutator(j1, j2, Word::new(terms))
                .expect("valid modified Q commutator");
        }
    }

    // 3. Power relations of N (shifted by +m)
    for i in 0..k {
        let shifted = shift_word(n.power(i), m);
        builder = builder
            .add_power(m + i, shifted)
            .expect("valid N power relation");
    }

    // 4. Commutator relations of N (shifted by +m)
    for i1 in 1..k {
        for i2 in 0..i1 {
            let shifted = shift_word(n.commutator(i1, i2), m);
            builder = builder
                .add_commutator(m + i1, m + i2, shifted)
                .expect("valid N commutator relation");
        }
    }

    // 5. Cross relations: all generators of N commute with Q ([g_{m+i}, g_j] = 1, default)

    let pres = builder.build();
    verify_consistency(&pres)?;
    Ok(pres)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zoo::{abelian, cyclic, dihedral};
    use crate::{chief_series, is_nilpotent, is_supersolvable, nilpotency_class};

    #[test]
    fn test_split_extension_s3() {
        let c3 = cyclic(3);
        let c2 = cyclic(2);

        // Action: q0 inverts n0 (n0 |-> n0^2)
        let action = vec![vec![Word::from_term(0, 2)]];
        let s3 = split_extension(&c3, &c2, &action).unwrap();

        assert_eq!(s3.order(), 6);
        assert_eq!(s3.num_gens(), 2);
        assert_eq!(verify_consistency(&s3), Ok(()));
        assert!(!is_abelian(&s3));
        assert!(is_supersolvable(&s3));
        assert!(!is_nilpotent(&s3));
    }

    #[test]
    fn test_split_extension_d8() {
        let c4 = cyclic(4);
        let c2 = cyclic(2);

        // Action: q0 inverts rotation generator n0 (n0 |-> n0^3 = n0 * n1)
        // c4 has generators n0 (order 2) and n1 (order 2) where n0^2 = n1.
        // Inverting n0 gives n0^3 = n0 n1.
        // Inverting n1 = n0^2 gives (n0^3)^2 = n0^6 = n0^2 = n1.
        let action = vec![vec![
            Word::from(vec![Term::new(0, 1), Term::new(1, 1)]),
            Word::from_term(1, 1),
        ]];
        let d8 = split_extension(&c4, &c2, &action).unwrap();

        assert_eq!(d8.order(), 8);
        assert_eq!(d8.num_gens(), 3);
        assert_eq!(verify_consistency(&d8), Ok(()));
        assert!(!is_abelian(&d8));
        assert!(is_nilpotent(&d8));
        assert_eq!(nilpotency_class(&d8), Some(2));
    }

    #[test]
    fn test_split_extension_a4() {
        let v4 = abelian(&[2, 2]);
        let c3 = cyclic(3);

        // Action of c3 (q0 of order 3) on V_4 = <n0, n1>:
        // q0 * n0 * q0^-1 = n1
        // q0 * n1 * q0^-1 = n0 * n1
        let action = vec![vec![
            Word::from_term(1, 1),
            Word::from(vec![Term::new(0, 1), Term::new(1, 1)]),
        ]];
        let a4 = split_extension(&v4, &c3, &action).unwrap();

        assert_eq!(a4.order(), 12);
        assert_eq!(a4.num_gens(), 3);
        assert_eq!(verify_consistency(&a4), Ok(()));
        assert!(!is_abelian(&a4));
        assert!(!is_nilpotent(&a4));
        assert!(!is_supersolvable(&a4));

        let chief = chief_series(&a4);
        let chief_orders: Vec<u128> = chief.iter().map(|h| h.order()).collect();
        assert_eq!(chief_orders, vec![12, 4, 1]);
    }

    #[test]
    fn test_split_extension_frobenius_f20() {
        let c5 = cyclic(5);
        let c4 = cyclic(4);

        // In F_20 ≅ C_5 ⋊ C_4, generator q0 (order 4, represented by [q0, q1] where q0^2 = q1)
        // acts on C_5 (n0) by multiplication by 2 (automorphism of order 4):
        // q0 * n0 * q0^-1 = n0^2
        // q1 * n0 * q1^-1 = (n0^2)^2 = n0^4
        let action = vec![vec![Word::from_term(0, 2)], vec![Word::from_term(0, 4)]];
        let f20 = split_extension(&c5, &c4, &action).unwrap();

        assert_eq!(f20.order(), 20);
        assert_eq!(f20.num_gens(), 3);
        assert_eq!(verify_consistency(&f20), Ok(()));
        assert!(!is_abelian(&f20));
        assert!(!is_nilpotent(&f20));
        assert!(is_supersolvable(&f20));
    }

    #[test]
    fn test_split_extension_invalid_action_fails_consistency() {
        let c5 = cyclic(5);
        let c2 = cyclic(2);

        // n0 |-> n0^2 has order 4 in Aut(C_5), which does not divide |C_2| = 2.
        // Therefore this action is not a valid homomorphism and fails confluence.
        let bad_action = vec![vec![Word::from_term(0, 2)]];
        let res = split_extension(&c5, &c2, &bad_action);
        assert!(matches!(res, Err(ExtensionError::Inconsistent(_))));
    }

    #[test]
    fn test_split_extension_direct_product() {
        let c2 = cyclic(2);
        let c3 = cyclic(3);

        // Trivial action yields direct product C_2 x C_3 ≅ C_6
        let c6 = split_extension(&c2, &c3, &[]).unwrap();
        assert_eq!(c6.order(), 6);
        assert_eq!(verify_consistency(&c6), Ok(()));
        assert!(is_abelian(&c6));
        assert!(is_nilpotent(&c6));
        assert_eq!(nilpotency_class(&c6), Some(1));
    }

    #[test]
    fn test_split_extension_non_abelian_kernel() {
        let d8 = dihedral(4);
        let c2 = cyclic(2);

        // Direct product D_8 x C_2
        let g = split_extension(&d8, &c2, &[]).unwrap();
        assert_eq!(g.order(), 16);
        assert_eq!(verify_consistency(&g), Ok(()));
        assert!(!is_abelian(&g));
        assert!(is_nilpotent(&g));
        assert_eq!(nilpotency_class(&g), Some(2));
    }

    #[test]
    fn test_split_extension_trivial_groups() {
        let c2 = cyclic(2);
        let trivial = Builder::new(vec![]).unwrap().build();

        // 1 ⋊ C_2 ≅ C_2
        let g1 = split_extension(&trivial, &c2, &[]).unwrap();
        assert_eq!(g1.order(), 2);
        assert_eq!(verify_consistency(&g1), Ok(()));

        // C_2 ⋊ 1 ≅ C_2
        let g2 = split_extension(&c2, &trivial, &[]).unwrap();
        assert_eq!(g2.order(), 2);
        assert_eq!(verify_consistency(&g2), Ok(()));

        // 1 ⋊ 1 ≅ 1
        let g3 = split_extension(&trivial, &trivial, &[]).unwrap();
        assert_eq!(g3.order(), 1);
        assert_eq!(verify_consistency(&g3), Ok(()));
    }

    #[test]
    fn test_central_extension_q8() {
        let c2 = cyclic(2);
        let v4 = abelian(&[2, 2]);

        let pow_tails = vec![Word::from_term(0, 1), Word::from_term(0, 1)];
        let comm_tails = vec![Word::from_term(0, 1)];

        let q8 = central_extension(&c2, &v4, &pow_tails, &comm_tails).unwrap();
        assert_eq!(q8.order(), 8);
        assert_eq!(verify_consistency(&q8), Ok(()));
        assert!(!is_abelian(&q8));
        assert!(is_nilpotent(&q8));
        assert_eq!(nilpotency_class(&q8), Some(2));
    }

    #[test]
    fn test_repeated_central_extensions_are_nilpotent() {
        // Repeated central extensions of cyclic groups of prime order always produce nilpotent groups.
        let c2 = cyclic(2);

        // Step 1: C_2 (order 2, class 1)
        let mut current_group = c2.clone();
        assert!(is_nilpotent(&current_group));
        assert_eq!(nilpotency_class(&current_group), Some(1));

        // Step 2: Central extension of C_2 by C_2 -> C_4 (order 4, class 1)
        let pow_tails_1 = vec![Word::from_term(0, 1)];
        let g2 = central_extension(&c2, &current_group, &pow_tails_1, &[]).unwrap();
        assert_eq!(g2.order(), 4);
        assert_eq!(verify_consistency(&g2), Ok(()));
        assert!(is_nilpotent(&g2));
        current_group = g2;

        // Step 3: Central extension of C_2 by C_4 -> order 8
        let pow_tails_2 = vec![Word::from_term(0, 1), Word::identity()];
        let g3 = central_extension(&c2, &current_group, &pow_tails_2, &[]).unwrap();
        assert_eq!(g3.order(), 8);
        assert_eq!(verify_consistency(&g3), Ok(()));
        assert!(is_nilpotent(&g3));
        current_group = g3;

        // Step 4: Central extension of C_2 by G_3 -> order 16
        let pow_tails_3 = vec![Word::identity(); current_group.num_gens()];
        let g4 = central_extension(&c2, &current_group, &pow_tails_3, &[]).unwrap();
        assert_eq!(g4.order(), 16);
        assert_eq!(verify_consistency(&g4), Ok(()));
        assert!(is_nilpotent(&g4));

        // Step 5: Test with p = 3 (Heisenberg / extraspecial group of order 27)
        let c3 = cyclic(3);
        let v9 = abelian(&[3, 3]);
        let comm_tails_p3 = vec![Word::from_term(0, 1)]; // [g1, g0] = g2 (central)
        let heisenberg = central_extension(&c3, &v9, &[], &comm_tails_p3).unwrap();
        assert_eq!(heisenberg.order(), 27);
        assert_eq!(verify_consistency(&heisenberg), Ok(()));
        assert!(!is_abelian(&heisenberg));
        assert!(is_nilpotent(&heisenberg));
        assert_eq!(nilpotency_class(&heisenberg), Some(2));
    }

    #[test]
    fn test_central_extension_non_abelian_kernel_rejected() {
        let d8 = dihedral(4);
        let c2 = cyclic(2);

        let res = central_extension(&d8, &c2, &[], &[]);
        assert_eq!(res, Err(ExtensionError::NonAbelianKernel));
    }

    #[test]
    fn test_extension_dimension_mismatch() {
        let c2 = cyclic(2);
        let c3 = cyclic(3);

        // Split extension dimension mismatch on Q
        let action_bad_q = vec![vec![Word::identity()], vec![Word::identity()]];
        assert_eq!(
            split_extension(&c2, &c3, &action_bad_q),
            Err(ExtensionError::DimensionMismatch {
                expected: 1,
                actual: 2,
            })
        );

        // Split extension dimension mismatch on N
        let action_bad_n = vec![vec![Word::identity(), Word::identity()]];
        assert_eq!(
            split_extension(&c2, &c3, &action_bad_n),
            Err(ExtensionError::DimensionMismatch {
                expected: 1,
                actual: 2,
            })
        );

        // Central extension dimension mismatch on powers
        assert_eq!(
            central_extension(&c2, &c3, &[Word::identity(), Word::identity()], &[]),
            Err(ExtensionError::DimensionMismatch {
                expected: 1,
                actual: 2,
            })
        );

        // Central extension dimension mismatch on commutators
        assert_eq!(
            central_extension(&c2, &c3, &[], &[Word::identity()]),
            Err(ExtensionError::DimensionMismatch {
                expected: 0,
                actual: 1,
            })
        );
    }

    #[test]
    fn test_extension_generator_out_of_bounds() {
        let c2 = cyclic(2);
        let c3 = cyclic(3);

        let bad_action = vec![vec![Word::from_term(5, 1)]];
        assert_eq!(
            split_extension(&c2, &c3, &bad_action),
            Err(ExtensionError::GeneratorOutOfBounds { index: 5, max: 0 })
        );

        let bad_tail = vec![Word::from_term(3, 1)];
        assert_eq!(
            central_extension(&c2, &c3, &bad_tail, &[]),
            Err(ExtensionError::GeneratorOutOfBounds { index: 3, max: 0 })
        );
    }

    #[test]
    fn test_error_display() {
        assert_eq!(
            ExtensionError::NonAbelianKernel.to_string(),
            "Central extension requires an abelian kernel group N"
        );
        assert_eq!(
            ExtensionError::DimensionMismatch {
                expected: 2,
                actual: 3
            }
            .to_string(),
            "Dimension mismatch: expected 2 elements, got 3"
        );
        assert_eq!(
            ExtensionError::GeneratorOutOfBounds { index: 4, max: 2 }.to_string(),
            "Generator index 4 is out of bounds (maximum valid index is 2)"
        );
    }
}
