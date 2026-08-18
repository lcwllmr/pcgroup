use crate::util::factorize;
use crate::{Builder, Presentation, Word};

/// Constructs a polycyclic (PC) presentation for the cyclic group `C_n` of order `n`.
///
/// If `n` factors into primes as `n = p_0 * p_1 * ... * p_{k-1}`, the presentation is given by:
/// `< g0, ..., g{k-1} | gi^pi = g{i+1} (0 <= i < k-1), g{k-1}^p{k-1} = 1, [gi, gj] = 1 (i > j) >`
///
/// For `n = 1`, returns the presentation of the trivial group with 0 generators.
///
/// # Examples
/// ```
/// use pcgroup::zoo::cyclic;
///
/// let c6 = cyclic(6);
/// assert_eq!(c6.order(), 6);
/// assert_eq!(c6.num_gens(), 2);
/// assert_eq!(c6.to_string(), "< g0, g1 | g0^2 = g1, g1^3 = 1, [g1, g0] = 1 >");
/// ```
///
/// # Panics
/// Panics if `n == 0`.
pub fn cyclic(n: u32) -> Presentation {
    assert!(
        n > 0,
        "Cyclic group order must be strictly positive (n >= 1)"
    );
    let factors = factorize(n);
    let mut builder =
        Builder::new(factors).expect("factorize produces valid prime relative orders");
    let k = builder.num_gens();
    for i in 0..k.saturating_sub(1) {
        builder = builder
            .add_power(i, Word::from_term(i + 1, 1))
            .expect("valid generator chain");
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
    fn test_cyclic_groups() {
        let orders = [1, 2, 3, 4, 5, 6, 8, 12, 16, 24, 60, 97, 360];
        for &n in &orders {
            let pres = cyclic(n);
            assert_eq!(
                pres.order(),
                n as u128,
                "Order mismatch for cyclic group C_{n}"
            );
            assert_eq!(
                verify_consistency(&pres),
                Ok(()),
                "Consistency check failed for cyclic group C_{n}"
            );

            // Cyclic groups are always abelian
            assert!(is_abelian(&pres), "Cyclic group C_{n} must be abelian");

            // Nilpotency class: 0 for trivial group C_1, 1 for C_n (n >= 2)
            let expected_class = if n == 1 { 0 } else { 1 };
            assert_eq!(
                nilpotency_class(&pres),
                Some(expected_class),
                "Nilpotency class mismatch for cyclic group C_{n}"
            );
            assert!(is_nilpotent(&pres), "Cyclic group C_{n} must be nilpotent");
            assert!(
                is_supersolvable(&pres),
                "Cyclic group C_{n} must be supersolvable"
            );

            // Irreducible representations check
            let irreps = irreducible_representations(&pres).expect("supersolvable");
            assert_eq!(irreps.len(), n as usize, "Irrep count mismatch for C_{n}");
            let sum_sq: usize = irreps.iter().map(|r| r.dim * r.dim).sum();
            assert_eq!(
                sum_sq, n as usize,
                "Sum of squared dimensions mismatch for C_{n}"
            );
            for rep in &irreps {
                assert_eq!(
                    rep.dim, 1,
                    "All irreps of cyclic group C_{n} must be 1-dimensional"
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "Cyclic group order must be strictly positive")]
    fn test_cyclic_zero_panics() {
        cyclic(0);
    }
}
