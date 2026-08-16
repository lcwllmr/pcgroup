//! Word collection engine for polycyclic (PC) presentations.
//!
//! The collection algorithm rewrites arbitrary generator words into canonical normal form
//! words using the power and commutator relations of a [`Presentation`].

use crate::word::Term;
use crate::{Presentation, Word};

/// An efficient stack-based word collector for Polycyclic Presentations.
///
/// # Algorithm
/// The collector maintains a collected prefix as an array of generator exponents `(e0, ..., e{n-1})`
/// where each `0 <= ei < pi`, and an uncollected suffix as a LIFO stack of [`Term`]s.
///
/// 1. If an incoming generator `gj` has `j >= max_active_gen(prefix)`, it is appended directly to `ej`,
///    triggering a power reduction if `ej >= pj`.
/// 2. If an incoming generator `gj` has `j < k` where `gk` is the largest active generator in the prefix,
///    one instance of `gk` is decremented from the prefix, and the commutation relation:
///    `gk * gj = gj * gk * [gk, gj]`
///    is applied by pushing the commutator tail `[gk, gj]`, `gk^1`, and `gj^1` onto the stack.
///
/// # Example: Collecting in S_3
/// ```
/// use pcgroup::{Builder, Collector, Word};
///
/// // Presentation of S_3: < g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >
/// let pres = Builder::new(vec![2, 3])
///     .unwrap()
///     .add_commutator(1, 0, Word::from_term(1, 1))
///     .unwrap()
///     .build();
///
/// // Collect the product g1 * g0 into normal form (g0 * g1^2)
/// let mut collector = Collector::new(&pres);
/// collector.collect(&Word::from_term(1, 1));
/// collector.collect(&Word::from_term(0, 1));
///
/// let result = collector.into_word();
/// assert_eq!(result.to_string(), "g0 g1^2");
/// ```
pub struct Collector<'a> {
    pres: &'a Presentation,

    /// The collected prefix represented as an exponent vector of length `n`.
    exponents: Vec<u32>,

    /// The uncollected suffix.
    /// Treated as a stack: the NEXT generator to process is at the end of the `Vec`.
    stack: Vec<Term>,
}

impl<'a> Collector<'a> {
    /// Initializes a new collector for the given presentation.
    pub fn new(pres: &'a Presentation) -> Self {
        Self {
            pres,
            exponents: vec![0; pres.num_gens()],
            // 64 is large enough to prevent reallocations for small/medium groups
            stack: Vec::with_capacity(64),
        }
    }

    /// Returns a slice of the collected prefix exponent vector.
    #[inline]
    pub fn exponents(&self) -> &[u32] {
        &self.exponents
    }

    /// Returns `true` if the collector currently represents the identity (empty prefix and stack).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.stack.is_empty() && self.exponents.iter().all(|&exp| exp == 0)
    }

    /// Clears the collected prefix and stack, resetting the collector for reuse.
    #[inline]
    pub fn reset(&mut self) {
        self.exponents.fill(0);
        self.stack.clear();
    }

    /// Finds the largest generator index currently active in the collected prefix.
    #[inline(always)]
    fn max_active_gen(&self) -> Option<usize> {
        // Reverse scan is fast for small n (generator count)
        self.exponents.iter().rposition(|&exp| exp > 0)
    }

    /// Collects a single generator power `gi^exponent` into the collected prefix.
    pub fn collect_generator(&mut self, gen_index: usize, exponent: u32) {
        if exponent == 0 {
            return;
        }
        self.collect_term(Term::new(gen_index, exponent));
    }

    /// Collects a single `Term` into the collected prefix.
    pub fn collect_term(&mut self, term: Term) {
        if term.exponent == 0 {
            return;
        }
        self.stack.push(term);
        self.process_stack();
    }

    /// Processes a normal form `Word` and incorporates it into the collected prefix.
    pub fn collect(&mut self, word: &Word) {
        if word.is_identity() {
            return;
        }

        // Fast-path: if collector prefix is empty, copy normal form word directly
        if self.is_empty() {
            for term in word.iter() {
                self.exponents[term.gen_index] = term.exponent;
            }
            return;
        }

        // Push the word onto the stack in reverse order so the first
        // generator of the word is processed first.
        for &term in word.iter().rev() {
            self.stack.push(term);
        }

        self.process_stack();
    }

    /// Core loop: processes all terms currently on the uncollected stack.
    fn process_stack(&mut self) {
        while let Some(term) = self.stack.pop() {
            let j = term.gen_index;
            let exp = term.exponent;

            if exp == 0 {
                continue;
            }

            if let Some(k) = self.max_active_gen()
                && k > j
            {
                // We must pass gj leftwards through gk.
                // Extract ONE gk from the prefix.
                self.exponents[k] -= 1;

                // If we had more than one gj, leave the rest on the stack for later
                if exp > 1 {
                    self.stack.push(Term::new(j, exp - 1));
                }

                // Relation: gk * gj = gj * gk * [gk, gj]
                let tail = self.pres.commutator(k, j);

                // 1. Push commutator tail in reverse order
                for &t in tail.iter().rev() {
                    self.stack.push(t);
                }

                // 2. Push gk^1
                self.stack.push(Term::new(k, 1));

                // 3. Push gj^1
                self.stack.push(Term::new(j, 1));

                continue;
            }

            // If we reach here, k <= j (or prefix is empty).
            // gj commutes with nothing in the prefix, so it safely appends.
            self.exponents[j] += exp;
            let p = self.pres.relative_order(j);

            // Check for power relations
            if self.exponents[j] >= p {
                let q = self.exponents[j] / p;
                self.exponents[j] %= p;

                let tail = self.pres.power(j);
                if !tail.is_empty() {
                    // Push `q` copies of the power tail, in reverse order
                    for _ in 0..q {
                        for &t in tail.iter().rev() {
                            self.stack.push(t);
                        }
                    }
                }
            }
        }
    }

    /// Converts the current exponent array state back into a collected normal form `Word`.
    pub fn into_word(self) -> Word {
        let active_count = self.exponents.iter().filter(|&&e| e > 0).count();
        let mut terms = Vec::with_capacity(active_count);

        for (g, &exp) in self.exponents.iter().enumerate() {
            if exp > 0 {
                terms.push(Term::new(g, exp));
            }
        }
        Word::new(terms)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Builder;

    #[test]
    fn test_collector_identity() {
        let pres = Builder::new(vec![2, 3]).unwrap().build();
        let mut collector = Collector::new(&pres);

        assert!(collector.is_empty());
        collector.collect(&Word::identity());
        assert!(collector.is_empty());

        let word = collector.into_word();
        assert!(word.is_identity());
    }

    #[test]
    fn test_collector_fast_path() {
        let pres = Builder::new(vec![2, 3, 5]).unwrap().build();
        let mut collector = Collector::new(&pres);

        let input = Word::from(vec![Term::new(0, 1), Term::new(2, 4)]);
        collector.collect(&input);
        assert_eq!(collector.exponents(), &[1, 0, 4]);

        let result = collector.into_word();
        assert_eq!(result, input);
    }

    #[test]
    fn test_collector_power_reductions() {
        let pres = Builder::new(vec![2, 3]).unwrap().build();
        let mut collector = Collector::new(&pres);

        // g0^4 in C_2 -> 1
        collector.collect_generator(0, 4);
        assert_eq!(collector.into_word(), Word::identity());

        // g1^5 in C_3 -> g1^2
        let mut collector2 = Collector::new(&pres);
        collector2.collect_generator(1, 5);
        assert_eq!(collector2.into_word(), Word::from_term(1, 2));
    }

    #[test]
    fn test_collector_power_cascade_c4() {
        // C_4: < g0, g1 | g0^2 = g1, g1^2 = 1, [g1, g0] = 1 >
        let pres = Builder::new(vec![2, 2])
            .unwrap()
            .add_power(0, Word::from_term(1, 1))
            .unwrap()
            .build();

        // g0 * g0 = g0^2 -> g1
        let mut col1 = Collector::new(&pres);
        col1.collect_generator(0, 1);
        col1.collect_generator(0, 1);
        assert_eq!(col1.into_word(), Word::from_term(1, 1));

        // g0^4 = 1
        let mut col2 = Collector::new(&pres);
        col2.collect_generator(0, 4);
        assert_eq!(col2.into_word(), Word::identity());

        // g0^3 = g0 * g1
        let mut col3 = Collector::new(&pres);
        col3.collect_generator(0, 3);
        assert_eq!(
            col3.into_word(),
            Word::from(vec![Term::new(0, 1), Term::new(1, 1)])
        );
    }

    #[test]
    fn test_collector_s3_commutator() {
        // S_3: < g0, g1 | g0^2 = 1, g1^3 = 1, [g1, g0] = g1 >
        let pres = Builder::new(vec![2, 3])
            .unwrap()
            .add_commutator(1, 0, Word::from_term(1, 1))
            .unwrap()
            .build();

        // g1 * g0 = g0 * g1^2
        let mut col1 = Collector::new(&pres);
        col1.collect(&Word::from_term(1, 1));
        col1.collect(&Word::from_term(0, 1));
        assert_eq!(
            col1.into_word(),
            Word::from(vec![Term::new(0, 1), Term::new(1, 2)])
        );

        // g0 * g1 * g0 = g1^2
        let mut col2 = Collector::new(&pres);
        col2.collect(&Word::from_term(0, 1));
        col2.collect(&Word::from_term(1, 1));
        col2.collect(&Word::from_term(0, 1));
        assert_eq!(col2.into_word(), Word::from_term(1, 2));

        // (g0 * g1)^2 = 1
        let mut col3 = Collector::new(&pres);
        col3.collect(&Word::from(vec![Term::new(0, 1), Term::new(1, 1)]));
        col3.collect(&Word::from(vec![Term::new(0, 1), Term::new(1, 1)]));
        assert_eq!(col3.into_word(), Word::identity());
    }

    #[test]
    fn test_collector_reset() {
        let pres = Builder::new(vec![2, 3]).unwrap().build();
        let mut collector = Collector::new(&pres);

        collector.collect_generator(0, 1);
        assert!(!collector.is_empty());

        collector.reset();
        assert!(collector.is_empty());
        assert_eq!(collector.exponents(), &[0, 0]);

        collector.collect_generator(1, 2);
        assert_eq!(collector.into_word(), Word::from_term(1, 2));
    }
}
