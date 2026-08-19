//! Finite field arithmetic for `F_q` where `q = p^r`.
//!
//! Elements of `F_q` are represented as polynomials in `F_p[x]` of degree `< r`
//! modulo a monic irreducible polynomial `f(x)` of degree `r`.

use crate::util::{factorize, mod_inverse};

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

/// Representation of a finite field `F_q` of order `q = p^r`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteField {
    /// Characteristic prime `p`.
    pub p: u32,
    /// Extension degree `r >= 1`.
    pub r: usize,
    /// Field order `q = p^r`.
    pub q: u64,
    /// Monic irreducible modulus polynomial of degree `r` over `F_p`.
    pub poly: Vec<u32>,
    /// Primitive element generating `F_q^\times`.
    pub prim: Vec<u32>,
}

impl FiniteField {
    /// Constructs field data for `F_{p^r}` by finding an irreducible polynomial and primitive root.
    pub fn new(p: u32, r: usize) -> Self {
        let q = (p as u64).pow(r as u32);
        if r == 1 {
            let order = p - 1;
            let mut prime_divs = factorize(order);
            prime_divs.dedup();
            let mut prim_elem = vec![1];
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
                    prim_elem = vec![g];
                    break;
                }
            }
            return Self {
                p,
                r,
                q,
                poly: vec![0, 1],
                prim: prim_elem,
            };
        }

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

        let mut prim = Vec::new();
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
                prim = elem;
                break;
            }
        }
        assert!(
            !prim.is_empty(),
            "primitive element in F_{}^{} not found",
            p,
            r
        );

        Self {
            p,
            r,
            q,
            poly,
            prim,
        }
    }

    /// Multiplies two field elements modulo the field polynomial.
    #[inline]
    pub fn mul(&self, a: &[u32], b: &[u32]) -> Vec<u32> {
        poly_rem(poly_mul(a, b, self.p), &self.poly, self.p)
    }

    /// Raises a field element to a power modulo the field polynomial.
    #[inline]
    pub fn pow(&self, base: &[u32], exp: u64) -> Vec<u32> {
        poly_powmod(base, exp, &self.poly, self.p)
    }

    /// Computes the multiplicative inverse of a non-zero element in `F_q^\times`.
    #[inline]
    pub fn inv(&self, a: &[u32]) -> Vec<u32> {
        if self.q <= 2 {
            return vec![1];
        }
        self.pow(a, self.q - 2)
    }

    /// Subtracts field elements: `a - b`.
    pub fn sub(&self, a: &[u32], b: &[u32]) -> Vec<u32> {
        let max_len = a.len().max(b.len());
        let mut res = vec![0; max_len];
        for (i, r) in res.iter_mut().enumerate() {
            let va = a.get(i).copied().unwrap_or(0);
            let vb = b.get(i).copied().unwrap_or(0);
            *r = (va + self.p - (vb % self.p)) % self.p;
        }
        trim(&mut res);
        res
    }

    /// Returns the `l`-th basis element `x^l` as a polynomial.
    pub fn basis_element(&self, l: usize) -> Vec<u32> {
        if l == 0 {
            vec![1]
        } else {
            let mut v = vec![0; l + 1];
            v[l] = 1;
            v
        }
    }

    /// Multiplies basis elements `x^l * x^v` and reduces modulo the field polynomial.
    pub fn basis_product(&self, l: usize, v: usize) -> Vec<u32> {
        let mut prod = vec![0; l + v + 1];
        prod[l + v] = 1;
        poly_rem(prod, &self.poly, self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_finite_field_basic() {
        let f4 = FiniteField::new(2, 2);
        assert_eq!(f4.q, 4);
        assert_eq!(f4.poly, vec![1, 1, 1]); // x^2 + x + 1

        let a = f4.basis_element(1); // x
        let a_sq = f4.mul(&a, &a); // x^2 = x + 1
        assert_eq!(a_sq, vec![1, 1]);
        assert_eq!(f4.inv(&a), a_sq);
    }
}
