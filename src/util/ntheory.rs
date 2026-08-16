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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_prime() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(5));
        assert!(!is_prime(6));
        assert!(is_prime(7));
        assert!(!is_prime(8));
        assert!(!is_prime(9));
        assert!(!is_prime(10));
        assert!(is_prime(11));
        assert!(is_prime(13));
        assert!(!is_prime(15));
        assert!(is_prime(17));
        assert!(is_prime(19));
        assert!(!is_prime(21));
        assert!(is_prime(23));
        assert!(!is_prime(25));
        assert!(is_prime(29));
        assert!(is_prime(31));
        assert!(!is_prime(35));
        assert!(is_prime(97));
        assert!(!is_prime(100));
    }
}
