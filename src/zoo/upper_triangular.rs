use crate::util::{FiniteField, factorize, is_prime};
use crate::word::Term;
use crate::zoo::unipotent::{populate_unipotent_commutators, unipotent_gen_index};
use crate::{Builder, Presentation, Word};

/// Constructs a PC presentation for the upper triangular Borel subgroup `B(n, q)` of `GL(n, q)`.
///
/// The group `B(n, q)` consists of all invertible `n x n` upper-triangular matrices with entries in the finite
/// field `F_q` (`q = p^r`, `p` prime, `r >= 1`).
///
/// The group has the semidirect product structure `B(n, q) = U(n, q) \rtimes T(n, q)` where `U(n, q)` is the [unipotent
/// subgroup](crate::zoo::unipotent) of unit upper-triangular matrices and `T(n, q) \cong (F_q^\times)^n` is the diagonal torus.
///
/// The order of `B(n, q)` is `(q - 1)^n * q^(n * (n - 1) / 2) = (p^r - 1)^n * p^(r * n * (n - 1) / 2)`.
///
/// Under the polycyclic presentation:
/// - Generators `g0, ..., g{n * m - 1}` represent the diagonal torus `T(n, q)` where `q - 1 = r_0 * ... * r_{m-1}`
///   is the prime factorization of `q - 1`.
/// - Generators `g_{n * m}, ..., g_{N-1}` represent the unipotent subgroup `U(n, q)` with relative order `p`.
///
/// Relations:
/// - Torus powers: `d_{i, s}^{r_s} = d_{i, s+1}` for `0 <= s < m - 1`, and `d_{i, m-1}^{r_{m-1}} = 1`.
/// - Torus-Torus commutators: `[d_A, d_B] = 1`.
/// - Unipotent powers: `u^p = 1`.
/// - Unipotent-Unipotent commutators: standard `U(n, q)` commutator relations.
/// - Torus-Unipotent commutators: `d_{a, s}` scales row `a` by `\omega_s^{-1} - 1`, and `d_{b, s}` scales column `b` by `\omega_s - 1`.
///
/// # Examples
/// ```
/// use pcgroup::zoo::upper_triangular;
///
/// // B(1, 3) is isomorphic to F_3^\times = C_2 of order 2
/// let b1_3 = upper_triangular(1, 3, 1);
/// assert_eq!(b1_3.order(), 2);
///
/// // B(2, 2) is isomorphic to U(2, 2) = C_2 of order 2
/// let b2_2 = upper_triangular(2, 2, 1);
/// assert_eq!(b2_2.order(), 2);
///
/// // B(2, 3) of order (2)^2 * 3 = 12
/// let b2_3 = upper_triangular(2, 3, 1);
/// assert_eq!(b2_3.order(), 12);
/// ```
///
/// # Panics
/// Panics if `n == 0`, if `p` is not prime, if `r == 0`, or if `p^r` overflows `u32`.
pub fn upper_triangular(n: usize, p: u32, r: u32) -> Presentation {
    assert!(n >= 1, "Dimension n must be strictly positive (n >= 1)");
    assert!(is_prime(p), "Field characteristic p must be prime (p >= 2)");
    assert!(r >= 1, "Field degree r must be strictly positive (r >= 1)");
    let q = p
        .checked_pow(r)
        .expect("Field size q = p^r exceeds u32 range");

    let r_usize = r as usize;
    let h_factors = factorize(q - 1);
    let m = h_factors.len();

    let num_t_gens = n * m;
    let num_pairs = n * (n - 1) / 2;
    let num_u_gens = num_pairs * r_usize;

    let mut rel_orders = Vec::with_capacity(num_t_gens + num_u_gens);
    for _ in 0..n {
        rel_orders.extend_from_slice(&h_factors);
    }
    rel_orders.extend(std::iter::repeat_n(p, num_u_gens));

    let mut builder = Builder::new(rel_orders).expect("valid prime relative orders");

    // 1. Torus power relations
    for i in 0..n {
        for s in 0..m.saturating_sub(1) {
            let cur = i * m + s;
            let next = cur + 1;
            builder = builder
                .add_power(cur, Word::from_term(next, 1))
                .expect("valid torus power relation");
        }
    }

    let field = FiniteField::new(p, r_usize);

    // 2. Unipotent commutators (shifted by num_t_gens)
    builder = populate_unipotent_commutators(builder, n, p, r_usize, num_t_gens, &field);

    // 3. Torus acting on Unipotent commutators
    if m > 0 && num_pairs > 0 {
        let mut alphas = Vec::with_capacity(m);
        let mut alphas_inv = Vec::with_capacity(m);
        let mut s_acc = 1u64;
        for &f in &h_factors {
            let elem = field.pow(&field.prim, s_acc);
            let elem_inv = field.inv(&elem);
            alphas.push(elem);
            alphas_inv.push(elem_inv);
            s_acc *= f as u64;
        }

        for x in 0..n {
            for s in 0..m {
                let t_idx = x * m + s;
                let alpha_s = &alphas[s];
                let alpha_s_inv = &alphas_inv[s];

                for a in 0..n {
                    for b in (a + 1)..n {
                        if x == a {
                            let diff = field.sub(alpha_s_inv, &[1]);
                            for l in 0..r_usize {
                                let x_l = field.basis_element(l);
                                let prod = field.mul(&diff, &x_l);
                                let mut terms = Vec::new();
                                for (w, &coeff) in prod.iter().enumerate() {
                                    let c = coeff % p;
                                    if c > 0 && w < r_usize {
                                        let target =
                                            num_t_gens + unipotent_gen_index(n, r_usize, a, b, w);
                                        terms.push(Term::new(target, c));
                                    }
                                }
                                if !terms.is_empty() {
                                    let u_idx =
                                        num_t_gens + unipotent_gen_index(n, r_usize, a, b, l);
                                    builder = builder
                                        .add_commutator(u_idx, t_idx, Word::new(terms))
                                        .expect("valid torus-unipotent commutator");
                                }
                            }
                        } else if x == b {
                            let diff = field.sub(alpha_s, &[1]);
                            for l in 0..r_usize {
                                let x_l = field.basis_element(l);
                                let prod = field.mul(&diff, &x_l);
                                let mut terms = Vec::new();
                                for (w, &coeff) in prod.iter().enumerate() {
                                    let c = coeff % p;
                                    if c > 0 && w < r_usize {
                                        let target =
                                            num_t_gens + unipotent_gen_index(n, r_usize, a, b, w);
                                        terms.push(Term::new(target, c));
                                    }
                                }
                                if !terms.is_empty() {
                                    let u_idx =
                                        num_t_gens + unipotent_gen_index(n, r_usize, a, b, l);
                                    builder = builder
                                        .add_commutator(u_idx, t_idx, Word::new(terms))
                                        .expect("valid torus-unipotent commutator");
                                }
                            }
                        }
                    }
                }
            }
        }
    }

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
    fn test_upper_triangular_groups() {
        let cases: &[(usize, u32, u32)] = &[
            (1, 2, 1),
            (1, 3, 1),
            (1, 2, 2),
            (1, 5, 1),
            (1, 7, 1),
            (1, 3, 2),
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
            (3, 5, 1),
            (4, 2, 1),
            (4, 3, 1),
        ];

        for &(n, p, r) in cases {
            let q = (p as u128).pow(r);
            let num_pairs = (n * (n - 1) / 2) as u32;
            let expected_order = (q - 1).pow(n as u32) * q.pow(num_pairs);
            let pres = upper_triangular(n, p, r);

            assert_eq!(
                pres.order(),
                expected_order,
                "Order mismatch for B({n}, {p}^{r}) = B({n}, {q})"
            );
            assert_eq!(
                verify_consistency(&pres),
                Ok(()),
                "Consistency check failed for B({n}, {p}^{r})"
            );

            // Abelian iff n == 1 or (n == 2 and q == 2)
            let is_ab = is_abelian(&pres);
            let expected_ab = n == 1 || (n == 2 && q == 2);
            assert_eq!(
                is_ab, expected_ab,
                "Abelian check mismatch for B({n}, {p}^{r})"
            );

            // Nilpotent iff n == 1 or q == 2
            let is_nil = is_nilpotent(&pres);
            let expected_nil = n == 1 || q == 2;
            assert_eq!(
                is_nil, expected_nil,
                "Nilpotency check mismatch for B({n}, {p}^{r})"
            );

            if expected_nil {
                let expected_class = if expected_order == 1 {
                    0
                } else if is_ab {
                    1
                } else {
                    n - 1
                };
                assert_eq!(
                    nilpotency_class(&pres),
                    Some(expected_class),
                    "Nilpotency class mismatch for B({n}, {p}^{r})"
                );
            } else {
                assert_eq!(
                    nilpotency_class(&pres),
                    None,
                    "Non-nilpotent B({n}, {p}^{r}) must return None for nilpotency_class"
                );
            }

            // Supersolvable iff n == 1 or r == 1 (prime field)
            let is_ss = is_supersolvable(&pres);
            let expected_ss = n == 1 || r == 1;
            assert_eq!(
                is_ss, expected_ss,
                "Supersolvability check mismatch for B({n}, {p}^{r}): n={n}, r={r}"
            );

            // Irrep test for supersolvable groups with moderate size
            if expected_ss && expected_order <= 300 {
                let irreps = irreducible_representations(&pres).expect("supersolvable");
                let sum_sq: usize = irreps.iter().map(|rep| rep.dim * rep.dim).sum();
                assert_eq!(
                    sum_sq, expected_order as usize,
                    "Sum of squared dimensions mismatch for B({n}, {p}^{r})"
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "Dimension n must be strictly positive")]
    fn test_upper_triangular_zero_n_panics() {
        upper_triangular(0, 2, 1);
    }

    #[test]
    #[should_panic(expected = "Field characteristic p must be prime")]
    fn test_upper_triangular_composite_p_panics() {
        upper_triangular(2, 4, 1);
    }

    #[test]
    #[should_panic(expected = "Field degree r must be strictly positive")]
    fn test_upper_triangular_zero_r_panics() {
        upper_triangular(2, 3, 0);
    }
}
