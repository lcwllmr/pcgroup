//! Exact arithmetic for finite-order monomial matrices over roots of unity.
//!
//! A monomial matrix has exactly one non-zero entry per row and per column. For an `e`-monomial
//! matrix, all non-zero entries are `e`-th roots of unity `exp(2*pi*i*vals[i]/e)`, represented
//! compactly by a permutation `perm` mapping row `i` to column `perm[i]` and an integer vector
//! `vals` of root-of-unity exponents in `Z/eZ`.

/// An `e`-monomial matrix whose non-zero entries are `e`-th roots of unity.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct MonomialMatrix {
    /// Permutation mapping row `i` to non-zero column `perm[i]`.
    pub perm: Vec<usize>,
    /// Root-of-unity phase exponents `vals[i]` representing entry `(i, perm[i]) = exp(2*pi*i*vals[i]/e)`.
    pub vals: Vec<u32>,
    /// Shared root-of-unity exponent `e`.
    pub e: u32,
}

impl MonomialMatrix {
    /// Constructs a new `MonomialMatrix` from a column permutation, phase values, and root exponent `e`.
    pub fn new(perm: Vec<usize>, vals: Vec<u32>, e: u32) -> Self {
        assert_eq!(
            perm.len(),
            vals.len(),
            "Permutation and values length mismatch"
        );
        let norm_vals = vals.into_iter().map(|v| v % e).collect();
        Self {
            perm,
            vals: norm_vals,
            e,
        }
    }

    /// Creates an identity `e`-monomial matrix of dimension `dim`.
    pub fn identity(dim: usize, e: u32) -> Self {
        Self {
            perm: (0..dim).collect(),
            vals: vec![0; dim],
            e,
        }
    }

    /// Creates a diagonal matrix with all diagonal entries equal to `exp(2*pi*i*val/e)`.
    pub fn diagonal(dim: usize, val: u32, e: u32) -> Self {
        Self {
            perm: (0..dim).collect(),
            vals: vec![val % e; dim],
            e,
        }
    }

    /// Creates an `e`-monomial matrix from a pure permutation matrix with all non-zero entries equal to `1`.
    pub fn from_perm(perm: &[usize], e: u32) -> Self {
        Self {
            perm: perm.to_vec(),
            vals: vec![0; perm.len()],
            e,
        }
    }

    /// Returns the dimension of the matrix.
    #[inline]
    pub fn dim(&self) -> usize {
        self.perm.len()
    }

    /// Returns the root-of-unity exponent `e`.
    #[inline]
    pub fn exponent(&self) -> u32 {
        self.e
    }

    /// Computes the product `self * other`.
    pub fn compose(&self, other: &Self) -> Self {
        debug_assert_eq!(self.dim(), other.dim());
        debug_assert_eq!(self.e, other.e);
        let dim = self.dim();
        let mut new_perm = Vec::with_capacity(dim);
        let mut new_vals = Vec::with_capacity(dim);
        for i in 0..dim {
            let j = self.perm[i];
            new_perm.push(other.perm[j]);
            new_vals.push((self.vals[i] + other.vals[j]) % self.e);
        }
        Self {
            perm: new_perm,
            vals: new_vals,
            e: self.e,
        }
    }

    /// Computes `self * other` storing the result directly into `out`.
    pub fn compose_into(&self, other: &Self, out: &mut Self) {
        debug_assert_eq!(self.dim(), other.dim());
        debug_assert_eq!(self.dim(), out.dim());
        debug_assert_eq!(self.e, other.e);
        debug_assert_eq!(self.e, out.e);
        let dim = self.dim();
        for i in 0..dim {
            let j = self.perm[i];
            out.perm[i] = other.perm[j];
            out.vals[i] = (self.vals[i] + other.vals[j]) % self.e;
        }
    }

    /// Computes the inverse matrix `self^-1`.
    pub fn invert(&self) -> Self {
        let dim = self.dim();
        let mut inv_perm = vec![0; dim];
        let mut inv_vals = vec![0; dim];
        for i in 0..dim {
            let j = self.perm[i];
            inv_perm[j] = i;
            inv_vals[j] = (self.e - (self.vals[i] % self.e)) % self.e;
        }
        Self {
            perm: inv_perm,
            vals: inv_vals,
            e: self.e,
        }
    }

    /// Writes the inverse permutation `i -> inv_perm[i]` into `out`.
    pub fn invert_perm_into(&self, out: &mut [usize]) {
        for (i, &j) in self.perm.iter().enumerate() {
            out[j] = i;
        }
    }

    /// Returns the inverse permutation vector `i -> inv_perm[i]`.
    pub fn inv_perm(&self) -> Vec<usize> {
        let mut inv = vec![0; self.dim()];
        self.invert_perm_into(&mut inv);
        inv
    }

    /// Computes the matrix power `self^p` for any integer `p` using binary exponentiation.
    pub fn pow(&self, mut p: i64) -> Self {
        let dim = self.dim();
        if p == 0 {
            return Self::identity(dim, self.e);
        }
        let base = if p < 0 {
            p = -p;
            self.invert()
        } else {
            self.clone()
        };

        let mut result = Self::identity(dim, self.e);
        let mut factor = base;
        let mut exp = p as u64;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.compose(&factor);
            }
            if exp > 1 {
                factor = factor.compose(&factor);
            }
            exp >>= 1;
        }
        result
    }

    /// Concatenates two monomial matrices into a block diagonal matrix `diag(self, other)`.
    pub fn direct_sum(&self, other: &Self) -> Self {
        debug_assert_eq!(self.e, other.e);
        let dim1 = self.dim();
        let dim2 = other.dim();
        let mut new_perm = Vec::with_capacity(dim1 + dim2);
        new_perm.extend_from_slice(&self.perm);
        for &p in &other.perm {
            new_perm.push(p + dim1);
        }
        let mut new_vals = Vec::with_capacity(dim1 + dim2);
        new_vals.extend_from_slice(&self.vals);
        new_vals.extend_from_slice(&other.vals);
        Self {
            perm: new_perm,
            vals: new_vals,
            e: self.e,
        }
    }

    /// Multiplies the matrix by a scalar root of unity `exp(2*pi*i*scalar/e)`.
    pub fn scalar_mul(&self, scalar: u32) -> Self {
        let norm_s = scalar % self.e;
        let vals = self.vals.iter().map(|&v| (v + norm_s) % self.e).collect();
        Self {
            perm: self.perm.clone(),
            vals,
            e: self.e,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_and_diagonal() {
        let id = MonomialMatrix::identity(4, 6);
        assert_eq!(id.dim(), 4);
        assert_eq!(id.exponent(), 6);
        assert_eq!(id.perm, vec![0, 1, 2, 3]);
        assert_eq!(id.vals, vec![0, 0, 0, 0]);

        let diag = MonomialMatrix::diagonal(3, 2, 6);
        assert_eq!(diag.dim(), 3);
        assert_eq!(diag.vals, vec![2, 2, 2]);
        assert_eq!(diag.pow(3), MonomialMatrix::identity(3, 6));
    }

    #[test]
    fn test_invert_and_compose() {
        let m = MonomialMatrix::new(vec![2, 0, 1], vec![1, 2, 3], 5);
        let inv = m.invert();
        assert_eq!(m.compose(&inv), MonomialMatrix::identity(3, 5));
        assert_eq!(inv.compose(&m), MonomialMatrix::identity(3, 5));
    }

    #[test]
    fn test_compose_into_and_inv_perm() {
        let m1 = MonomialMatrix::new(vec![2, 0, 1], vec![1, 2, 3], 5);
        let m2 = MonomialMatrix::new(vec![1, 2, 0], vec![4, 3, 2], 5);
        let mut out = MonomialMatrix::identity(3, 5);
        m1.compose_into(&m2, &mut out);
        assert_eq!(out, m1.compose(&m2));
        assert_eq!(m1.inv_perm(), vec![1, 2, 0]);
    }

    #[test]
    fn test_powers() {
        let m = MonomialMatrix::new(vec![1, 2, 0], vec![0, 0, 1], 4);
        assert_eq!(m.pow(0), MonomialMatrix::identity(3, 4));
        assert_eq!(m.pow(1), m);
        assert_eq!(m.pow(3), m.compose(&m).compose(&m));
        assert_eq!(m.pow(-1), m.invert());
        assert_eq!(m.pow(-2), m.invert().compose(&m.invert()));
    }

    #[test]
    fn test_direct_sum_and_scalar_mul() {
        let m1 = MonomialMatrix::new(vec![1, 0], vec![1, 2], 4);
        let m2 = MonomialMatrix::new(vec![0], vec![3], 4);
        let sum = m1.direct_sum(&m2);
        assert_eq!(sum.dim(), 3);
        assert_eq!(sum.perm, vec![1, 0, 2]);
        assert_eq!(sum.vals, vec![1, 2, 3]);

        let scaled = m1.scalar_mul(2);
        assert_eq!(scaled.vals, vec![3, 0]);
    }

    #[test]
    fn test_from_perm() {
        let perm = [2, 0, 1];
        let m = MonomialMatrix::from_perm(&perm, 7);
        assert_eq!(m.dim(), 3);
        assert_eq!(m.vals, vec![0, 0, 0]);
        assert_eq!(m.perm, vec![2, 0, 1]);
    }
}
