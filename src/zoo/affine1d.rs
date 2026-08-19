use crate::util::{FiniteField, factorize, is_prime};
use crate::word::Term;
use crate::{Builder, Presentation, Word};

/// Constructs a polycyclic (PC) presentation for the 1-dimensional affine general linear group `AGL(1, q)`.
///
/// The group is defined as the affine transformations `x |-> a * x + b` of the finite field `F_q`
/// where `q = p^r` (`p` is prime, `r >= 1`), with `a in F_q^\times` and `b in F_q`.
///
/// The group structure is the semidirect product `AGL(1, q) = N \rtimes H \cong (F_q, +) \rtimes F_q^\times \cong (C_p)^r \rtimes C_{q-1}`,
/// of order `q * (q - 1) = p^r * (p^r - 1)`.
///
/// Under the polycyclic presentation:
/// - Generators `g0, ..., g{m-1}` represent the multiplicative subgroup `H \cong C_{q-1}` with prime relative orders
///   `f_0, ..., f_{m-1}` where `q - 1 = f_0 * ... * f_{m-1}`.
/// - Generators `g_m, ..., g_{m+r-1}` represent the translation subgroup `N \cong (C_p)^r` with prime relative orders `p`.
///
/// Relations:
/// - `g_j^{f_j} = g_{j+1}` for `0 <= j < m - 1`, and `g_{m-1}^{f_{m-1}} = 1`
/// - `g_{m+i}^p = 1` for `0 <= i < r`
/// - `[g_i, g_j] = 1` for `0 <= j < i < m` and for `m <= j < i < m + r`
/// - `[g_{m+i}, g_j] = (\alpha_j - 1) * x^i` expressed in the polynomial basis of `F_q`
///
/// # Examples
/// ```
/// use pcgroup::zoo::affine1d;
///
/// // AGL(1, 2) = C_2 of order 2
/// let agl2 = affine1d(2, 1);
/// assert_eq!(agl2.order(), 2);
///
/// // AGL(1, 3) = S_3 of order 6
/// let s3 = affine1d(3, 1);
/// assert_eq!(s3.order(), 6);
/// assert_eq!(s3.to_string(), "< g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >");
///
/// // AGL(1, 4) = A_4 of order 12
/// let a4 = affine1d(2, 2);
/// assert_eq!(a4.order(), 12);
/// ```
///
/// # Panics
/// Panics if `p` is not prime, if `r == 0`, or if `p^r` overflows `u32`.
pub fn affine1d(p: u32, r: u32) -> Presentation {
    assert!(is_prime(p), "Field characteristic p must be prime (p >= 2)");
    assert!(r >= 1, "Exponent r must be strictly positive (r >= 1)");
    let q = p
        .checked_pow(r)
        .expect("Group order q = p^r exceeds u32 range");

    if q == 2 {
        return Builder::new(vec![2]).expect("2 is prime").build();
    }

    let h_factors = factorize(q - 1);
    let m = h_factors.len();
    let r_usize = r as usize;

    let mut rel_orders = h_factors.clone();
    rel_orders.extend(std::iter::repeat_n(p, r_usize));

    let mut builder = Builder::new(rel_orders).expect("valid prime relative orders");

    for j in 0..m.saturating_sub(1) {
        builder = builder
            .add_power(j, Word::from_term(j + 1, 1))
            .expect("valid H power relation");
    }

    let field = FiniteField::new(p, r_usize);

    let mut alphas = Vec::with_capacity(m);
    let mut s = 1u64;
    for &f in &h_factors {
        alphas.push(field.pow(&field.prim, s));
        s *= f as u64;
    }

    for i in 0..r_usize {
        let x_i = field.basis_element(i);

        for (j, alpha_j) in alphas.iter().enumerate() {
            let alpha_minus_1 = field.sub(alpha_j, &[1]);
            let prod = field.mul(&alpha_minus_1, &x_i);
            let mut terms = Vec::new();
            for (k, &coeff) in prod.iter().enumerate() {
                if coeff > 0 && k < r_usize {
                    terms.push(Term::new(m + k, coeff));
                }
            }
            builder = builder
                .add_commutator(m + i, j, Word::new(terms))
                .expect("valid commutator relation");
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
    fn test_affine1d_groups() {
        let cases: &[(u32, u32)] = &[
            (2, 1),
            (2, 2),
            (2, 3),
            (2, 4),
            (3, 1),
            (3, 2),
            (3, 3),
            (5, 1),
            (5, 2),
            (7, 1),
            (11, 1),
            (17, 1),
        ];

        for &(p, r) in cases {
            let q = (p as u128).pow(r);
            let expected_order = q * (q - 1);
            let pres = affine1d(p, r);

            assert_eq!(
                pres.order(),
                expected_order,
                "Order mismatch for AGL(1, {p}^{r}) = AGL(1, {q})"
            );
            assert_eq!(
                verify_consistency(&pres),
                Ok(()),
                "Consistency check failed for AGL(1, {p}^{r})"
            );

            let is_ab = is_abelian(&pres);
            assert_eq!(is_ab, q == 2, "Abelian check mismatch for AGL(1, {p}^{r})");

            let is_nil = is_nilpotent(&pres);
            assert_eq!(
                is_nil,
                q == 2,
                "Nilpotency check mismatch for AGL(1, {p}^{r})"
            );

            let expected_class = if q == 2 { Some(1) } else { None };
            assert_eq!(
                nilpotency_class(&pres),
                expected_class,
                "Nilpotency class mismatch for AGL(1, {p}^{r})"
            );

            // AGL(1, p^r) is supersolvable if and only if r == 1
            let is_ss = is_supersolvable(&pres);
            assert_eq!(
                is_ss,
                r == 1,
                "Supersolvability mismatch for AGL(1, {p}^{r}): r={r}"
            );

            // Irreducible representations check for supersolvable cases (r = 1)
            if r == 1 {
                let irreps = irreducible_representations(&pres).expect("supersolvable");
                let count_1d = irreps.iter().filter(|rep| rep.dim == 1).count();
                let expected_1d = if p == 2 { 2 } else { (p - 1) as usize };
                assert_eq!(
                    count_1d, expected_1d,
                    "1D irrep count mismatch for AGL(1, {p})"
                );

                if p > 2 {
                    let count_high_d = irreps
                        .iter()
                        .filter(|rep| rep.dim == (p - 1) as usize)
                        .count();
                    assert_eq!(
                        count_high_d, 1,
                        "High-dimensional irrep count mismatch for AGL(1, {p})"
                    );
                }

                assert_eq!(
                    irreps.len(),
                    if p == 2 { 2 } else { p as usize },
                    "Total irrep count mismatch for AGL(1, {p})"
                );
                let sum_sq: usize = irreps.iter().map(|rep| rep.dim * rep.dim).sum();
                assert_eq!(
                    sum_sq,
                    (p * (p - 1)) as usize,
                    "Sum of squared dimensions mismatch for AGL(1, {p})"
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "Field characteristic p must be prime")]
    fn test_affine1d_composite_p_panics() {
        affine1d(4, 1);
    }

    #[test]
    #[should_panic(expected = "Field characteristic p must be prime")]
    fn test_affine1d_zero_p_panics() {
        affine1d(0, 1);
    }

    #[test]
    #[should_panic(expected = "Exponent r must be strictly positive")]
    fn test_affine1d_zero_r_panics() {
        affine1d(3, 0);
    }
}
