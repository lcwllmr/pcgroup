/// Returns `true` if `n` is a prime number (`n >= 2`).
#[inline]
pub fn is_prime(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n.is_multiple_of(2) || n.is_multiple_of(3) {
        return false;
    }
    let mut d = 5u32;
    while (d as u64) * (d as u64) <= n as u64 {
        if n.is_multiple_of(d) || n.is_multiple_of(d + 2) {
            return false;
        }
        d += 6;
    }
    true
}

/// Returns the prime factors of `n` in non-decreasing order with multiplicity.
///
/// For `n <= 1`, an empty vector is returned.
pub fn factorize(mut n: u32) -> Vec<u32> {
    let mut factors = Vec::new();
    if n <= 1 {
        return factors;
    }
    while n.is_multiple_of(2) {
        factors.push(2);
        n /= 2;
    }
    while n.is_multiple_of(3) {
        factors.push(3);
        n /= 3;
    }
    let mut d = 5u32;
    while (d as u64) * (d as u64) <= n as u64 {
        while n.is_multiple_of(d) {
            factors.push(d);
            n /= d;
        }
        while n.is_multiple_of(d + 2) {
            factors.push(d + 2);
            n /= d + 2;
        }
        d += 6;
    }
    if n > 1 {
        factors.push(n);
    }
    factors
}

/// Computes the modular multiplicative inverse of `a` modulo `m` such that `(a * x) % m == 1`.
///
/// Returns `None` if `m <= 1` or if `gcd(a, m) != 1`.
pub fn mod_inverse(a: u32, m: u32) -> Option<u32> {
    if m <= 1 {
        return None;
    }
    let mut t = 0i64;
    let mut newt = 1i64;
    let mut r = m as i64;
    let mut newr = (a % m) as i64;

    while newr != 0 {
        let q = r / newr;
        let next_t = t - q * newt;
        t = newt;
        newt = next_t;
        let next_r = r - q * newr;
        r = newr;
        newr = next_r;
    }

    if r > 1 {
        return None;
    }
    if t < 0 {
        t += m as i64;
    }
    Some(t as u32)
}

/// Computes the greatest common divisor of `a` and `b`.
#[inline]
pub fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Computes the least common multiple of `a` and `b`.
#[inline]
pub fn lcm(a: u32, b: u32) -> u32 {
    if a == 0 || b == 0 {
        0
    } else {
        ((a as u64 / gcd(a, b) as u64) * b as u64) as u32
    }
}

/// Solves `p * c == s (mod e)` for all `c in 0..e`.
///
/// In the context of roots of unity, this computes all `p`-th roots of `exp(2*pi*i*s/e)`
/// in the cyclic group of `e`-th roots of unity.
pub fn extract_roots(s: u32, p: u32, e: u32) -> Vec<u32> {
    if e == 0 {
        return Vec::new();
    }
    let s_norm = s % e;
    let g = gcd(p, e);
    if !s_norm.is_multiple_of(g) {
        return Vec::new();
    }
    let p_prime = p / g;
    let s_prime = s_norm / g;
    let e_prime = e / g;

    if e_prime == 1 {
        let mut roots = Vec::with_capacity(g as usize);
        for k in 0..g {
            roots.push(k);
        }
        return roots;
    }

    let inv = mod_inverse(p_prime, e_prime).unwrap_or(0);
    let c0 = ((s_prime as u64 * inv as u64) % e_prime as u64) as u32;

    let mut roots = Vec::with_capacity(g as usize);
    for k in 0..g {
        roots.push(c0 + k * e_prime);
    }
    roots
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_prime() {
        let primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 97];
        let non_primes = [0, 1, 4, 6, 8, 9, 10, 15, 21, 25, 35, 100];
        for &p in &primes {
            assert!(is_prime(p), "{p} should be prime");
        }
        for &np in &non_primes {
            assert!(!is_prime(np), "{np} should not be prime");
        }
    }

    #[test]
    fn test_factorize() {
        let cases: &[(u32, &[u32])] = &[
            (0, &[]),
            (1, &[]),
            (2, &[2]),
            (3, &[3]),
            (4, &[2, 2]),
            (12, &[2, 2, 3]),
            (60, &[2, 2, 3, 5]),
            (97, &[97]),
            (360, &[2, 2, 2, 3, 3, 5]),
        ];
        for &(n, expected) in cases {
            assert_eq!(factorize(n), expected, "factorization mismatch for {n}");
        }
    }

    #[test]
    fn test_mod_inverse() {
        assert_eq!(mod_inverse(3, 7), Some(5)); // 3 * 5 = 15 = 1 mod 7
        assert_eq!(mod_inverse(2, 5), Some(3)); // 2 * 3 = 6 = 1 mod 5
        assert_eq!(mod_inverse(1, 13), Some(1));
        assert_eq!(mod_inverse(0, 5), None);
        assert_eq!(mod_inverse(2, 4), None);
        assert_eq!(mod_inverse(5, 1), None);
    }

    #[test]
    fn test_extract_roots() {
        assert_eq!(extract_roots(0, 2, 4), vec![0, 2]);
        assert_eq!(extract_roots(2, 2, 4), vec![1, 3]);
        assert_eq!(extract_roots(0, 3, 6), vec![0, 2, 4]);
        assert_eq!(extract_roots(1, 2, 4), vec![]);
    }
}
