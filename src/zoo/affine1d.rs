use crate::util::{factorize, is_prime, mod_inverse};
use crate::word::Term;
use crate::{Builder, Presentation, Word};

fn trim(poly: &mut Vec<u32>) {
    while poly.len() > 1 && *poly.last().unwrap() == 0 {
        poly.pop();
    }
}

fn poly_mul(a: &[u32], b: &[u32], p: u32) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return vec![0];
    }
    let mut res = vec![0; a.len() + b.len() - 1];
    for (i, &va) in a.iter().enumerate() {
        for (j, &vb) in b.iter().enumerate() {
            res[i + j] = ((res[i + j] as u64 + va as u64 * vb as u64) % p as u64) as u32;
        }
    }
    trim(&mut res);
    res
}

fn poly_rem(mut a: Vec<u32>, b: &[u32], p: u32) -> Vec<u32> {
    if a.len() < b.len() {
        trim(&mut a);
        return a;
    }
    let inv_lead = mod_inverse(*b.last().unwrap() % p, p).unwrap_or(1);
    while a.len() >= b.len() {
        let lead_a = *a.last().unwrap();
        if lead_a != 0 {
            let factor = ((lead_a as u64 * inv_lead as u64) % p as u64) as u32;
            let shift = a.len() - b.len();
            for (i, &vb) in b.iter().enumerate() {
                let sub = ((vb as u64 * factor as u64) % p as u64) as u32;
                a[shift + i] = (a[shift + i] + p - sub) % p;
            }
        }
        a.pop();
    }
    trim(&mut a);
    if a.is_empty() { vec![0] } else { a }
}

fn poly_gcd(mut a: Vec<u32>, mut b: Vec<u32>, p: u32) -> Vec<u32> {
    while !(b.len() == 1 && b[0] == 0) && !b.is_empty() {
        let r = poly_rem(a, &b, p);
        a = b;
        b = r;
    }
    a
}

fn poly_powmod(base: &[u32], mut exp: u64, modulus: &[u32], p: u32) -> Vec<u32> {
    let mut res = vec![1];
    let mut cur = poly_rem(base.to_vec(), modulus, p);
    while exp > 0 {
        if exp & 1 == 1 {
            res = poly_rem(poly_mul(&res, &cur, p), modulus, p);
        }
        cur = poly_rem(poly_mul(&cur, &cur, p), modulus, p);
        exp >>= 1;
    }
    res
}

fn find_field(p: u32, r: usize) -> (Vec<u32>, Vec<u32>) {
    if r == 1 {
        let order = p - 1;
        let mut prime_divs = factorize(order);
        prime_divs.dedup();
        for g in 1..p {
            let mut is_prim = true;
            for &d in &prime_divs {
                let mut cur = 1u64;
                let mut base = g as u64;
                let mut exp = (order / d) as u64;
                while exp > 0 {
                    if exp & 1 == 1 {
                        cur = (cur * base) % (p as u64);
                    }
                    base = (base * base) % (p as u64);
                    exp >>= 1;
                }
                if cur == 1 {
                    is_prim = false;
                    break;
                }
            }
            if is_prim {
                return (vec![0, 1], vec![g]);
            }
        }
        return (vec![0, 1], vec![1]);
    }

    let q = (p as u64).pow(r as u32);
    let mut r_prime_divs = factorize(r as u32);
    r_prime_divs.dedup();

    let mut poly = vec![0; r + 1];
    poly[r] = 1;
    let mut found_poly = false;
    for code in 1..q {
        let mut val = code;
        for coeff in poly.iter_mut().take(r) {
            *coeff = (val % (p as u64)) as u32;
            val /= p as u64;
        }
        if poly[0] == 0 {
            continue;
        }

        let x = vec![0, 1];
        if poly_powmod(&x, q, &poly, p) != x {
            continue;
        }

        let mut is_irred = true;
        for &d in &r_prime_divs {
            let exp_d = (p as u64).pow((r as u32) / d);
            let mut diff = poly_powmod(&x, exp_d, &poly, p);
            if diff.len() < 2 {
                diff.resize(2, 0);
            }
            diff[1] = (diff[1] + p - 1) % p;
            trim(&mut diff);
            if poly_gcd(poly.clone(), diff, p).len() > 1 {
                is_irred = false;
                break;
            }
        }
        if is_irred {
            found_poly = true;
            break;
        }
    }
    assert!(
        found_poly,
        "monic irreducible polynomial of degree {r} over F_{p} not found"
    );

    let order = q - 1;
    let mut q_prime_divs = factorize(order as u32);
    q_prime_divs.dedup();

    for code in 1..q {
        let mut elem = vec![0; r];
        let mut val = code;
        for coeff in elem.iter_mut().take(r) {
            *coeff = (val % (p as u64)) as u32;
            val /= p as u64;
        }
        trim(&mut elem);
        if elem.len() == 1 && elem[0] == 0 {
            continue;
        }

        let mut is_prim = true;
        for &d in &q_prime_divs {
            let exp = order / (d as u64);
            if poly_powmod(&elem, exp, &poly, p) == vec![1] {
                is_prim = false;
                break;
            }
        }
        if is_prim {
            return (poly, elem);
        }
    }
    panic!("primitive element in F_{}^{} not found", p, r);
}

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

    let (poly, prim) = find_field(p, r_usize);

    let mut alphas = Vec::with_capacity(m);
    let mut s = 1u64;
    for &f in &h_factors {
        alphas.push(poly_powmod(&prim, s, &poly, p));
        s *= f as u64;
    }

    for i in 0..r_usize {
        let x_i = if i == 0 {
            vec![1]
        } else {
            let mut v = vec![0; i + 1];
            v[i] = 1;
            v
        };

        for (j, alpha_j) in alphas.iter().enumerate() {
            let mut alpha_minus_1 = alpha_j.clone();
            if alpha_minus_1.is_empty() {
                alpha_minus_1 = vec![0];
            }
            alpha_minus_1[0] = (alpha_minus_1[0] + p - 1) % p;
            trim(&mut alpha_minus_1);

            let prod = poly_rem(poly_mul(&alpha_minus_1, &x_i, p), &poly, p);
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
    use crate::{is_abelian, is_nilpotent, is_supersolvable, nilpotency_class, verify_consistency};

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
