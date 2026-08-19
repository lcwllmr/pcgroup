use crate::util::{FiniteField, is_prime};
use crate::word::Term;
use crate::{Builder, Presentation, Word};

/// Computes the unipotent pair index for `(i, j)` where `0 <= i < j < n`.
#[inline]
pub(crate) fn unipotent_pair_index(n: usize, i: usize, j: usize) -> usize {
    let h = j - i;
    ((h - 1) * n - (h - 1) * h / 2) + i
}

/// Computes the unipotent generator index for `(i, j, l)`.
#[inline]
pub(crate) fn unipotent_gen_index(n: usize, r: usize, i: usize, j: usize, l: usize) -> usize {
    unipotent_pair_index(n, i, j) * r + l
}

/// Populates unipotent commutator relations into the presentation builder.
pub(crate) fn populate_unipotent_commutators(
    mut builder: Builder,
    n: usize,
    p: u32,
    r: usize,
    offset: usize,
    field: &FiniteField,
) -> Builder {
    for i in 0..n {
        for j in (i + 1)..n {
            for y in (j + 1)..n {
                for l in 0..r {
                    for v in 0..r {
                        let i1 = offset + unipotent_gen_index(n, r, i, j, l);
                        let i2 = offset + unipotent_gen_index(n, r, j, y, v);
                        let prod = field.basis_product(l, v);

                        if i2 > i1 {
                            let mut terms = Vec::new();
                            for (w, &coeff) in prod.iter().enumerate() {
                                if coeff > 0 && w < r {
                                    let neg_c = (p - (coeff % p)) % p;
                                    if neg_c > 0 {
                                        let target = offset + unipotent_gen_index(n, r, i, y, w);
                                        terms.push(Term::new(target, neg_c));
                                    }
                                }
                            }
                            if !terms.is_empty() {
                                builder = builder
                                    .add_commutator(i2, i1, Word::new(terms))
                                    .expect("valid unipotent commutator");
                            }
                        } else {
                            let mut terms = Vec::new();
                            for (w, &coeff) in prod.iter().enumerate() {
                                let c = coeff % p;
                                if c > 0 && w < r {
                                    let target = offset + unipotent_gen_index(n, r, i, y, w);
                                    terms.push(Term::new(target, c));
                                }
                            }
                            if !terms.is_empty() {
                                builder = builder
                                    .add_commutator(i1, i2, Word::new(terms))
                                    .expect("valid unipotent commutator");
                            }
                        }
                    }
                }
            }
        }
    }
    builder
}

/// Constructs a PC presentation for the unipotent subgroup `U(n, q)` of `GL(n, q)`.
///
/// The group `U(n, q)` consists of all `n x n` unit upper-triangular matrices with entries in the finite
/// field `F_q` (`q = p^r`, `p` prime, `r >= 1`).
///
/// The order of `U(n, q)` is `q^(n * (n - 1) / 2) = p^(r * n * (n - 1) / 2)`.
///
/// Under the polycyclic presentation:
/// - Each generator corresponds to an elementary unipotent matrix `I + x^l * E_{i, j}` where
///   `0 <= i < j < n` and `0 <= l < r`.
/// - Generators are ordered primarily by ascending height `h = j - i`, secondarily by row `i`,
///   and tertiarily by basis exponent `l`.
/// - All generators have prime relative order `p` (`u^p = 1`).
///
/// # Examples
/// ```
/// use pcgroup::zoo::unipotent;
///
/// // U(1, q) is trivial of order 1
/// let u1 = unipotent(1, 3, 1);
/// assert_eq!(u1.order(), 1);
///
/// // U(2, 4) is elementary abelian (C_2)^2 of order 4
/// let u2_4 = unipotent(2, 2, 2);
/// assert_eq!(u2_4.order(), 4);
///
/// // U(3, 2) is the Heisenberg group of order 8 (isomorphic to D_8)
/// let u3_2 = unipotent(3, 2, 1);
/// assert_eq!(u3_2.order(), 8);
/// assert_eq!(u3_2.to_string(), "< g0, g1, g2 | g0^2 = 1, g1^2 = 1, g2^2 = 1, [g1, g0] = g2, [g2, g0] = 1, [g2, g1] = 1 >");
/// ```
///
/// # Panics
/// Panics if `n == 0`, if `p` is not prime, if `r == 0`, or if `p^r` overflows `u32`.
pub fn unipotent(n: usize, p: u32, r: u32) -> Presentation {
    assert!(n >= 1, "Dimension n must be strictly positive (n >= 1)");
    assert!(is_prime(p), "Field characteristic p must be prime (p >= 2)");
    assert!(r >= 1, "Field degree r must be strictly positive (r >= 1)");
    assert!(
        p.checked_pow(r).is_some(),
        "Field size q = p^r exceeds u32 range"
    );

    if n == 1 {
        return Builder::new(Vec::new())
            .expect("empty relative orders")
            .build();
    }

    let r_usize = r as usize;
    let num_pairs = n * (n - 1) / 2;
    let num_gens = num_pairs * r_usize;
    let rel_orders = vec![p; num_gens];

    let builder = Builder::new(rel_orders).expect("valid prime relative orders");
    let field = FiniteField::new(p, r_usize);

    let builder = populate_unipotent_commutators(builder, n, p, r_usize, 0, &field);
    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        irreducible_representations, is_abelian, is_nilpotent, is_supersolvable, nilpotency_class,
        verify_consistency,
    };

    #[test]
    fn test_unipotent_groups() {
        let cases: &[(usize, u32, u32)] = &[
            (1, 2, 1),
            (1, 3, 1),
            (1, 2, 2),
            (1, 5, 1),
            (2, 2, 1),
            (2, 2, 2),
            (2, 2, 3),
            (2, 3, 1),
            (2, 3, 2),
            (2, 5, 1),
            (2, 7, 1),
            (3, 2, 1),
            (3, 2, 2),
            (3, 3, 1),
            (3, 3, 2),
            (3, 5, 1),
            (3, 7, 1),
            (4, 2, 1),
            (4, 2, 2),
            (4, 3, 1),
            (5, 2, 1),
        ];

        for &(n, p, r) in cases {
            let q = (p as u128).pow(r);
            let num_pairs = (n * (n - 1) / 2) as u32;
            let expected_order = q.pow(num_pairs);
            let pres = unipotent(n, p, r);

            assert_eq!(
                pres.order(),
                expected_order,
                "Order mismatch for U({n}, {p}^{r}) = U({n}, {q})"
            );
            assert_eq!(
                verify_consistency(&pres),
                Ok(()),
                "Consistency check failed for U({n}, {p}^{r})"
            );

            // Abelian iff n <= 2
            let is_ab = is_abelian(&pres);
            assert_eq!(is_ab, n <= 2, "Abelian check mismatch for U({n}, {p}^{r})");

            // Always nilpotent (p-group)
            assert!(
                is_nilpotent(&pres),
                "Unipotent group U({n}, {p}^{r}) must be nilpotent"
            );

            // Nilpotency class: 0 for n=1, 1 for n=2, n-1 for n >= 3
            let expected_class = if n == 1 {
                0
            } else if n == 2 {
                1
            } else {
                n - 1
            };
            assert_eq!(
                nilpotency_class(&pres),
                Some(expected_class),
                "Nilpotency class mismatch for U({n}, {p}^{r})"
            );

            // Always supersolvable (finite p-group)
            assert!(
                is_supersolvable(&pres),
                "Unipotent group U({n}, {p}^{r}) must be supersolvable"
            );

            // Irrep test for moderately sized groups
            if expected_order <= 300 {
                let irreps = irreducible_representations(&pres).expect("supersolvable");
                let sum_sq: usize = irreps.iter().map(|rep| rep.dim * rep.dim).sum();
                assert_eq!(
                    sum_sq, expected_order as usize,
                    "Sum of squared dimensions mismatch for U({n}, {p}^{r})"
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "Dimension n must be strictly positive")]
    fn test_unipotent_zero_n_panics() {
        unipotent(0, 2, 1);
    }

    #[test]
    #[should_panic(expected = "Field characteristic p must be prime")]
    fn test_unipotent_composite_p_panics() {
        unipotent(2, 4, 1);
    }

    #[test]
    #[should_panic(expected = "Field degree r must be strictly positive")]
    fn test_unipotent_zero_r_panics() {
        unipotent(2, 3, 0);
    }
}
