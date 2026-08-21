use crate::word::Term;
use crate::zoo::{affine1d, dihedral};
use crate::{Builder, Presentation, Word};

/// Constructs a polycyclic (PC) presentation for the symmetric group `S_3` of order 6.
///
/// The symmetric group `S_3` consists of all 6 permutations of 3 elements.
/// It is isomorphic to the dihedral group `D_6` and to the 1-dimensional affine group `AGL(1, 3)`.
///
/// Under the polycyclic presentation:
/// - Generator `g0` represents the transposition `(1, 2)` (relative order 2).
/// - Generator `g1` represents the 3-cycle `(1, 2, 3)` (relative order 3).
///
/// Relations:
/// - `g0^2 = 1`
/// - `g1^3 = 1`
/// - `[g1, g0] = g1`
///
/// # Examples
/// ```
/// use pcgroup::zoo::s3;
///
/// let g = s3();
/// assert_eq!(g.order(), 6);
/// assert_eq!(g.num_gens(), 2);
/// assert_eq!(g.to_string(), "< g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >");
/// ```
pub fn s3() -> Presentation {
    dihedral(3)
}

/// Constructs a polycyclic (PC) presentation for the alternating group `A_4` of order 12.
///
/// The alternating group `A_4` consists of all 12 even permutations of 4 elements.
/// It is isomorphic to the rotational symmetry group of a regular tetrahedron and to
/// the 1-dimensional affine group `AGL(1, 4) ≅ V_4 ⋊ C_3`.
///
/// Under the polycyclic presentation:
/// - Generator `g0` represents a 3-cycle `(1, 2, 3)` (relative order 3).
/// - Generator `g1` represents the double transposition `(1, 2)(3, 4)` (relative order 2).
/// - Generator `g2` represents the double transposition `(1, 3)(2, 4)` (relative order 2).
///
/// Relations:
/// - `g0^3 = 1`, `g1^2 = 1`, `g2^2 = 1`
/// - `[g1, g0] = g1 g2`
/// - `[g2, g0] = g1`
/// - `[g2, g1] = 1`
///
/// # Examples
/// ```
/// use pcgroup::zoo::a4;
///
/// let g = a4();
/// assert_eq!(g.order(), 12);
/// assert_eq!(g.num_gens(), 3);
/// assert_eq!(
///     g.to_string(),
///     "< g0, g1, g2 | g0^3 = 1, g1^2 = 1, g2^2 = 1, [g1, g0] = g1 g2, [g2, g0] = g1, [g2, g1] = 1 >"
/// );
/// ```
pub fn a4() -> Presentation {
    affine1d(2, 2)
}

/// Constructs a polycyclic (PC) presentation for the symmetric group `S_4` of order 24.
///
/// The symmetric group `S_4` consists of all 24 permutations of 4 elements.
/// It is isomorphic to the full symmetry group of a regular octahedron / cube and
/// has the subnormal series `S_4 ▷ A_4 ▷ V_4 ▷ C_2 ▷ 1`.
///
/// Under the polycyclic presentation:
/// - Generator `g0` represents the transposition `(1, 2)` (relative order 2).
/// - Generator `g1` represents the 3-cycle `(1, 2, 3)` (relative order 3).
/// - Generator `g2` represents the double transposition `(1, 2)(3, 4)` (relative order 2).
/// - Generator `g3` represents the double transposition `(1, 3)(2, 4)` (relative order 2).
///
/// Relations:
/// - `g0^2 = 1`, `g1^3 = 1`, `g2^2 = 1`, `g3^2 = 1`
/// - `[g1, g0] = g1`
/// - `[g2, g0] = 1`
/// - `[g2, g1] = g2 g3`
/// - `[g3, g0] = g2`
/// - `[g3, g1] = g2`
/// - `[g3, g2] = 1`
///
/// # Examples
/// ```
/// use pcgroup::zoo::s4;
///
/// let g = s4();
/// assert_eq!(g.order(), 24);
/// assert_eq!(g.num_gens(), 4);
/// assert_eq!(
///     g.to_string(),
///     "< g0, g1, g2, g3 | g0^2 = 1, g1^3 = 1, g2^2 = 1, g3^2 = 1, [g1, g0] = g1, [g2, g0] = 1, [g2, g1] = g2 g3, [g3, g0] = g2, [g3, g1] = g2, [g3, g2] = 1 >"
/// );
/// ```
pub fn s4() -> Presentation {
    Builder::new(vec![2, 3, 2, 2])
        .expect("valid relative orders")
        .add_commutator(1, 0, Word::from_term(1, 1))
        .expect("valid relation")
        .add_commutator(3, 0, Word::from_term(2, 1))
        .expect("valid relation")
        .add_commutator(2, 1, Word::new(vec![Term::new(2, 1), Term::new(3, 1)]))
        .expect("valid relation")
        .add_commutator(3, 1, Word::from_term(2, 1))
        .expect("valid relation")
        .build()
}

/// Constructs a polycyclic (PC) presentation for the special linear group `SL(2, 3)` of order 24.
///
/// The group `SL(2, 3)` consists of all 2x2 matrices with entries in `F_3` and determinant 1.
/// It is isomorphic to the binary tetrahedral group `2T`, the double cover of `A_4 ≅ PSL(2, 3)`.
///
/// It has the semidirect product structure `SL(2, 3) ≅ Q_8 ⋊ C_3` where `C_3` cyclically permutes
/// the quaternion elements `j ↦ i ↦ k ↦ j`.
///
/// Under the polycyclic presentation:
/// - Generator `g0` represents an element of order 3 in `SL(2, 3)` (relative order 3).
/// - Generator `g1` represents quaternion generator `j` (relative order 2, `g1^2 = g3`).
/// - Generator `g2` represents quaternion generator `i` (relative order 2, `g2^2 = g3`).
/// - Generator `g3` represents the central matrix `-I` (relative order 2).
///
/// Relations:
/// - `g0^3 = 1`, `g1^2 = g3`, `g2^2 = g3`, `g3^2 = 1`
/// - `[g1, g0] = g1 g2 g3`
/// - `[g2, g0] = g1`
/// - `[g2, g1] = g3`
/// - `[g3, g0] = 1`
/// - `[g3, g1] = 1`
/// - `[g3, g2] = 1`
///
/// # Examples
/// ```
/// use pcgroup::zoo::sl2_3;
///
/// let g = sl2_3();
/// assert_eq!(g.order(), 24);
/// assert_eq!(g.num_gens(), 4);
/// assert_eq!(
///     g.to_string(),
///     "< g0, g1, g2, g3 | g0^3 = 1, g1^2 = g3, g2^2 = g3, g3^2 = 1, [g1, g0] = g1 g2 g3, [g2, g0] = g1, [g2, g1] = g3, [g3, g0] = 1, [g3, g1] = 1, [g3, g2] = 1 >"
/// );
/// ```
pub fn sl2_3() -> Presentation {
    Builder::new(vec![3, 2, 2, 2])
        .expect("valid relative orders")
        .add_power(1, Word::from_term(3, 1))
        .expect("valid relation")
        .add_power(2, Word::from_term(3, 1))
        .expect("valid relation")
        .add_commutator(
            1,
            0,
            Word::new(vec![Term::new(1, 1), Term::new(2, 1), Term::new(3, 1)]),
        )
        .expect("valid relation")
        .add_commutator(2, 0, Word::from_term(1, 1))
        .expect("valid relation")
        .add_commutator(2, 1, Word::from_term(3, 1))
        .expect("valid relation")
        .build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        GeneratingSequence, chief_series, commutator_subgroup, irreducible_representations,
        is_abelian, is_nilpotent, is_supersolvable, verify_consistency,
    };

    #[test]
    fn test_outliers_consistency_and_orders() {
        for (pres, expected_order, expected_gens) in
            [(s3(), 6, 2), (a4(), 12, 3), (s4(), 24, 4), (sl2_3(), 24, 4)]
        {
            assert_eq!(pres.order(), expected_order);
            assert_eq!(pres.num_gens(), expected_gens);
            assert_eq!(verify_consistency(&pres), Ok(()));
            assert!(!is_abelian(&pres));
            assert!(!is_nilpotent(&pres));
        }
    }

    #[test]
    fn test_supersolvability_and_chief_series() {
        // S_3 is supersolvable
        assert!(is_supersolvable(&s3()));
        let cs_s3: Vec<u128> = chief_series(&s3()).iter().map(|h| h.order()).collect();
        assert_eq!(cs_s3, vec![6, 3, 1]);

        // A_4 is solvable but not supersolvable (chief factor V_4 of order 4)
        assert!(!is_supersolvable(&a4()));
        let cs_a4: Vec<u128> = chief_series(&a4()).iter().map(|h| h.order()).collect();
        assert_eq!(cs_a4, vec![12, 4, 1]);

        // S_4 is solvable but not supersolvable (chief factor V_4 of order 4)
        assert!(!is_supersolvable(&s4()));
        let cs_s4: Vec<u128> = chief_series(&s4()).iter().map(|h| h.order()).collect();
        assert_eq!(cs_s4, vec![24, 12, 4, 1]);

        // SL(2, 3) is solvable but not supersolvable
        assert!(!is_supersolvable(&sl2_3()));
        let cs_sl2_3: Vec<u128> = chief_series(&sl2_3()).iter().map(|h| h.order()).collect();
        assert_eq!(cs_sl2_3, vec![24, 8, 2, 1]);
    }

    #[test]
    fn test_derived_subgroups() {
        // [S_3, S_3] = A_3 (order 3)
        let s3_pres = s3();
        let full_s3 = GeneratingSequence::full_group(&s3_pres);
        let s3_gens = full_s3.elements();
        assert_eq!(commutator_subgroup(&full_s3, &full_s3, s3_gens).order(), 3);

        // [A_4, A_4] = V_4 (order 4)
        let a4_pres = a4();
        let full_a4 = GeneratingSequence::full_group(&a4_pres);
        let a4_gens = full_a4.elements();
        assert_eq!(commutator_subgroup(&full_a4, &full_a4, a4_gens).order(), 4);

        // [S_4, S_4] = A_4 (order 12)
        let s4_pres = s4();
        let full_s4 = GeneratingSequence::full_group(&s4_pres);
        let s4_gens = full_s4.elements();
        assert_eq!(commutator_subgroup(&full_s4, &full_s4, s4_gens).order(), 12);

        // [SL(2, 3), SL(2, 3)] = Q_8 (order 8)
        let sl_pres = sl2_3();
        let full_sl = GeneratingSequence::full_group(&sl_pres);
        let sl_gens = full_sl.elements();
        assert_eq!(commutator_subgroup(&full_sl, &full_sl, sl_gens).order(), 8);
    }

    #[test]
    fn test_s3_irreps() {
        let irreps = irreducible_representations(&s3()).expect("s3 is supersolvable");
        let sum_sq: usize = irreps.iter().map(|r| r.dim * r.dim).sum();
        assert_eq!(sum_sq, 6);
        assert_eq!(irreps.iter().filter(|r| r.dim == 1).count(), 2);
        assert_eq!(irreps.iter().filter(|r| r.dim == 2).count(), 1);
    }
}
