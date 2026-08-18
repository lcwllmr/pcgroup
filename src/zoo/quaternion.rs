use crate::util::factorize;
use crate::word::Term;
use crate::zoo::dihedral::inv_sq_word;
use crate::{Builder, Presentation, Word};

/// Constructs a polycyclic (PC) presentation for the generalized quaternion group `Q_{4n}` of order `4n` (`n >= 2`).
///
/// The group is defined as:
/// `< x, y | x^{2n} = 1, y^2 = x^n, y^-1 * x * y = x^-1 >`
///
/// For `n = 2`, this constructs the standard quaternion group `Q_8` of order 8.
///
/// Under the polycyclic presentation, generator `g0` represents `y` (relative order 2),
/// and generators `g1, ..., gk` represent the cyclic subgroup `< x > = C_{2n}`
/// with prime relative orders `p_1, ..., p_k` (`2n = p_1 * ... * p_k`, where `p_1 = 2`).
///
/// Relations:
/// - `g0^2 = x^n` expressed in the polycyclic basis of `< x >`
/// - `gi^pi = g{i+1}` for `1 <= i < k`, and `gk^pk = 1`
/// - `[gi, gj] = 1` for all `1 <= j < i <= k`
/// - `[gi, g0] = gi^-2` expressed in normal form for `1 <= i <= k`
///
/// # Examples
/// ```
/// use pcgroup::zoo::quaternion;
///
/// // Quaternion group Q_8 of order 8 (n = 2)
/// let q8 = quaternion(2);
/// assert_eq!(q8.order(), 8);
/// assert_eq!(q8.num_gens(), 3);
/// assert_eq!(
///     q8.to_string(),
///     "< g0, g1, g2 | g0^2 = g2, g1^2 = g2, g2^2 = 1, [g1, g0] = g2, [g2, g0] = 1, [g2, g1] = 1 >"
/// );
///
/// // Dicyclic group / generalized quaternion Q_12 of order 12 (n = 3)
/// let q12 = quaternion(3);
/// assert_eq!(q12.order(), 12);
/// assert_eq!(q12.num_gens(), 3);
/// ```
///
/// # Panics
/// Panics if `n < 2`.
pub fn quaternion(n: u32) -> Presentation {
    assert!(
        n >= 2,
        "Quaternion group order parameter n must be at least 2 (order 4n >= 8)"
    );

    let mut rel_orders = vec![2];
    rel_orders.extend(factorize(2 * n));

    let mut builder = Builder::new(rel_orders.clone()).expect("valid prime relative orders");
    let k = builder.num_gens() - 1;

    // Power relation for g0: g0^2 = x^n in the basis of <x>
    let mut temp = n;
    let mut y_sq_terms = Vec::new();
    for (m, &p_m) in rel_orders.iter().enumerate().skip(1) {
        let digit = temp % p_m;
        if digit > 0 {
            y_sq_terms.push(Term::new(m, digit));
        }
        temp /= p_m;
    }
    builder = builder
        .add_power(0, Word::new(y_sq_terms))
        .expect("valid power relation for g0");

    // Power relations in the cyclic subgroup <x>
    for i in 1..k {
        builder = builder
            .add_power(i, Word::from_term(i + 1, 1))
            .expect("valid cyclic subgroup power relation");
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
    use crate::{is_abelian, is_nilpotent, is_supersolvable, nilpotency_class, verify_consistency};

    #[test]
    fn test_quaternion_groups() {
        let orders = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 30];
        for &n in &orders {
            let pres = quaternion(n);
            assert_eq!(
                pres.order(),
                (4 * n) as u128,
                "Order mismatch for quaternion group Q_{{{}}}",
                4 * n
            );
            assert_eq!(
                verify_consistency(&pres),
                Ok(()),
                "Consistency check failed for quaternion group Q_{{{}}}",
                4 * n
            );

            // Quaternion / dicyclic groups Q_{4n} (n >= 2) are never abelian
            assert!(
                !is_abelian(&pres),
                "Quaternion group Q_{{{}}} must not be abelian",
                4 * n
            );

            // Q_{4n} is nilpotent iff n is a power of 2
            let is_pwr2 = n.is_power_of_two();
            assert_eq!(
                is_nilpotent(&pres),
                is_pwr2,
                "Nilpotency check mismatch for quaternion group Q_{{{}}}",
                4 * n
            );

            if is_pwr2 {
                // Nilpotency class: k + 1 for Q_{2^(k+2)} where n = 2^k (e.g. Q_8 -> 2, Q_16 -> 3, Q_32 -> 4)
                let expected_class = (n.trailing_zeros() + 1) as usize;
                assert_eq!(
                    nilpotency_class(&pres),
                    Some(expected_class),
                    "Nilpotency class mismatch for quaternion group Q_{{{}}}",
                    4 * n
                );
            } else {
                assert_eq!(
                    nilpotency_class(&pres),
                    None,
                    "Non-2-group quaternion Q_{{{}}} must not be nilpotent",
                    4 * n
                );
            }

            // All generalized quaternion / dicyclic groups Q_{4n} are supersolvable
            assert!(
                is_supersolvable(&pres),
                "Quaternion group Q_{{{}}} must be supersolvable",
                4 * n
            );
        }
    }

    #[test]
    #[should_panic(expected = "Quaternion group order parameter n must be at least 2")]
    fn test_quaternion_underflow_panics() {
        quaternion(1);
    }
}
