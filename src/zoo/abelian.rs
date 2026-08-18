use crate::util::factorize;
use crate::{Builder, Presentation, Word};

/// Constructs a polycyclic (PC) presentation for a finite abelian group `C_{d_0} x ... x C_{d_{m-1}}`.
///
/// Each non-trivial direct cyclic factor `C_{d_r}` is decomposed into prime power cycles
/// with sequential chained power relations, while all commutator relations are trivial (`[gi, gj] = 1`).
///
/// # Examples
/// ```
/// use pcgroup::zoo::abelian;
///
/// // Klein four-group V_4 = C_2 x C_2
/// let v4 = abelian(&[2, 2]);
/// assert_eq!(v4.order(), 4);
/// assert_eq!(v4.num_gens(), 2);
/// assert_eq!(v4.to_string(), "< g0, g1 | g0^2 = 1, g1^2 = 1, [g1, g0] = 1 >");
///
/// // C_4 x C_6 of order 24
/// let g = abelian(&[4, 6]);
/// assert_eq!(g.order(), 24);
/// assert_eq!(g.num_gens(), 4);
/// ```
///
/// # Panics
/// Panics if any element of `orders` is zero.
pub fn abelian(orders: &[u32]) -> Presentation {
    let mut all_factors = Vec::new();
    let mut factor_ranges = Vec::new();

    for &order in orders {
        assert!(
            order > 0,
            "Direct factor orders must be strictly positive (d >= 1)"
        );
        let factors = factorize(order);
        if !factors.is_empty() {
            let start = all_factors.len();
            let count = factors.len();
            all_factors.extend(factors);
            factor_ranges.push(start..start + count);
        }
    }

    let mut builder =
        Builder::new(all_factors).expect("factorize produces valid prime relative orders");
    for range in factor_ranges {
        if range.len() > 1 {
            for i in range.start..range.end - 1 {
                builder = builder
                    .add_power(i, Word::from_term(i + 1, 1))
                    .expect("valid power relation within factor");
            }
        }
    }
    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{is_abelian, is_nilpotent, is_supersolvable, nilpotency_class, verify_consistency};

    #[test]
    fn test_abelian_groups() {
        let cases: &[(&[u32], u128)] = &[
            (&[], 1),
            (&[1, 1], 1),
            (&[2, 2], 4),
            (&[2, 3], 6),
            (&[4, 6], 24),
            (&[2, 2, 2], 8),
            (&[3, 3, 3], 27),
            (&[12, 18, 30], 6480),
        ];

        for &(factors, expected_order) in cases {
            let pres = abelian(factors);
            assert_eq!(
                pres.order(),
                expected_order,
                "Order mismatch for abelian group with factors {factors:?}"
            );
            assert_eq!(
                verify_consistency(&pres),
                Ok(()),
                "Consistency check failed for abelian group with factors {factors:?}"
            );

            // Direct products of cyclic groups are always abelian
            assert!(
                is_abelian(&pres),
                "Abelian group with factors {factors:?} must be abelian"
            );

            // Nilpotency class: 0 for trivial group, 1 for non-trivial abelian groups
            let expected_class = if expected_order == 1 { 0 } else { 1 };
            assert_eq!(
                nilpotency_class(&pres),
                Some(expected_class),
                "Nilpotency class mismatch for abelian group with factors {factors:?}"
            );
            assert!(
                is_nilpotent(&pres),
                "Abelian group with factors {factors:?} must be nilpotent"
            );
            assert!(
                is_supersolvable(&pres),
                "Abelian group with factors {factors:?} must be supersolvable"
            );
        }
    }

    #[test]
    #[should_panic(expected = "Direct factor orders must be strictly positive")]
    fn test_abelian_zero_panics() {
        abelian(&[2, 0, 3]);
    }
}
