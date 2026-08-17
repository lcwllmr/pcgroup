use crate::util::factorize;
use crate::word::Term;
use crate::{Builder, Presentation, Word};

/// Builds the normal form word for `gi^-2` in the cyclic subgroup `< g1, ..., gk >`.
pub(crate) fn inv_sq_word(i: usize, rel_orders: &[u32]) -> Word {
    let p_i = rel_orders[i];
    let mut terms = Vec::new();
    if p_i > 2 {
        terms.push(Term::new(i, p_i - 2));
    }
    for (m, &p_m) in rel_orders.iter().enumerate().skip(i + 1) {
        terms.push(Term::new(m, p_m - 1));
    }
    Word::new(terms)
}

/// Constructs a polycyclic (PC) presentation for the dihedral group `D_{2n}` of order `2n` (`n >= 1`).
///
/// The group is defined as the symmetries of a regular `n`-gon:
/// `< s, r | s^2 = 1, r^n = 1, s * r * s = r^-1 >`
///
/// Under the polycyclic presentation, generator `g0` represents the reflection `s` (relative order 2),
/// and generators `g1, ..., gk` represent the rotation subgroup `< r > = C_n`
/// with prime relative orders `p_1, ..., p_k` (`n = p_1 * ... * p_k`).
///
/// Relations:
/// - `g0^2 = 1`
/// - `gi^pi = g{i+1}` for `1 <= i < k`, and `gk^pk = 1`
/// - `[gi, gj] = 1` for all `1 <= j < i <= k`
/// - `[gi, g0] = gi^-2` expressed in normal form for `1 <= i <= k`
///
/// # Examples
/// ```
/// use pcgroup::zoo::dihedral;
///
/// // S_3 = D_6 of order 6
/// let s3 = dihedral(3);
/// assert_eq!(s3.order(), 6);
/// assert_eq!(s3.num_gens(), 2);
/// assert_eq!(s3.to_string(), "< g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >");
///
/// // D_8 of order 8
/// let d8 = dihedral(4);
/// assert_eq!(d8.order(), 8);
/// assert_eq!(d8.num_gens(), 3);
/// assert_eq!(
///     d8.to_string(),
///     "< g0, g1, g2 | g0^2 = 1, g1^2 = g2, g2^2 = 1, [g1, g0] = g2, [g2, g0] = 1, [g2, g1] = 1 >"
/// );
/// ```
///
/// # Panics
/// Panics if `n == 0`.
pub fn dihedral(n: u32) -> Presentation {
    assert!(
        n > 0,
        "Dihedral group parameter n must be strictly positive (n >= 1)"
    );
    if n == 1 {
        // D_2 = C_2
        return Builder::new(vec![2]).expect("2 is prime").build();
    }

    let mut rel_orders = vec![2];
    rel_orders.extend(factorize(n));

    let mut builder = Builder::new(rel_orders.clone()).expect("valid prime relative orders");
    let k = builder.num_gens() - 1;

    // Power relations in the rotation subgroup
    for i in 1..k {
        builder = builder
            .add_power(i, Word::from_term(i + 1, 1))
            .expect("valid rotation power relation");
    }

    // Commutator relations [gi, g0] = gi^-2
    for i in 1..=k {
        let comm_word = inv_sq_word(i, &rel_orders);
        builder = builder
            .add_commutator(i, 0, comm_word)
            .expect("valid commutator relation");
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Element, GeneratingSequence, is_abelian, is_nilpotent, nilpotency_class, verify_consistency,
    };

    #[test]
    fn test_dihedral_groups() {
        let orders = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 30, 60];
        for &n in &orders {
            let pres = dihedral(n);
            assert_eq!(
                pres.order(),
                (2 * n) as u128,
                "Order mismatch for dihedral group D_{{{}}}",
                2 * n
            );
            assert_eq!(
                verify_consistency(&pres),
                Ok(()),
                "Consistency check failed for dihedral group D_{{{}}}",
                2 * n
            );

            let full = GeneratingSequence::full_group(&pres);
            let gens: Vec<_> = (0..pres.num_gens())
                .map(|i| Element::from_generator(i, 1, &pres))
                .collect();

            // Dihedral D_{2n} is abelian iff n = 1 (D_2 = C_2) or n = 2 (D_4 = V_4)
            let is_ab = is_abelian(&full, &gens);
            assert_eq!(
                is_ab,
                n <= 2,
                "Abelian check mismatch for dihedral group D_{{{}}}",
                2 * n
            );

            // D_{2n} is nilpotent iff n is a power of 2
            let is_pwr2 = n.is_power_of_two();
            assert_eq!(
                is_nilpotent(&full, &gens),
                is_pwr2,
                "Nilpotency check mismatch for dihedral group D_{{{}}}",
                2 * n
            );

            if is_pwr2 {
                // Nilpotency class: 1 for D_2 and D_4; k = log2(n) for D_{2^(k+1)} (n = 2^k >= 4)
                let expected_class = if n <= 2 {
                    1
                } else {
                    n.trailing_zeros() as usize
                };
                assert_eq!(
                    nilpotency_class(&full, &gens),
                    Some(expected_class),
                    "Nilpotency class mismatch for dihedral group D_{{{}}}",
                    2 * n
                );
            } else {
                assert_eq!(
                    nilpotency_class(&full, &gens),
                    None,
                    "Non-2-group dihedral D_{{{}}} must not be nilpotent",
                    2 * n
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "Dihedral group parameter n must be strictly positive")]
    fn test_dihedral_zero_panics() {
        dihedral(0);
    }
}
