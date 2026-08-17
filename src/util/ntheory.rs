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
}
