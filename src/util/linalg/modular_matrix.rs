//! Finite-field linear algebra and module operations over `F_p`.
//!
//! This module provides matrix operations over finite fields `F_p`, row reduction,
//! and basic representation theory algorithms (a lite Spin Algorithm / MeatAxe) to find
//! submodule chains for chief series computation.

#![allow(clippy::needless_range_loop)]
#![allow(unused_assignments)]
#![allow(clippy::manual_is_multiple_of)]

use crate::util::mod_inverse;

/// A matrix over the finite field `F_p`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModularMatrix {
    pub data: Vec<Vec<u32>>,
    pub p: u32,
    pub rows: usize,
    pub cols: usize,
}

impl ModularMatrix {
    /// Creates a new matrix with all zeros.
    pub fn zeros(rows: usize, cols: usize, p: u32) -> Self {
        Self {
            data: vec![vec![0; cols]; rows],
            p,
            rows,
            cols,
        }
    }

    /// Creates an `n x n` identity matrix over `F_p`.
    pub fn identity(n: usize, p: u32) -> Self {
        let mut data = vec![vec![0; n]; n];
        for i in 0..n {
            data[i][i] = 1 % p;
        }
        Self {
            data,
            p,
            rows: n,
            cols: n,
        }
    }

    /// Creates a matrix from column vectors.
    pub fn from_columns(cols: &[Vec<u32>], p: u32) -> Self {
        let num_cols = cols.len();
        let num_rows = if num_cols > 0 { cols[0].len() } else { 0 };
        let mut data = vec![vec![0u32; num_cols]; num_rows];
        for (c, col) in cols.iter().enumerate() {
            for (r, &val) in col.iter().enumerate() {
                data[r][c] = val % p;
            }
        }
        Self {
            data,
            p,
            rows: num_rows,
            cols: num_cols,
        }
    }

    /// Generates the standard basis vectors for `F_p^d`.
    pub fn standard_basis(d: usize) -> Vec<Vec<u32>> {
        let mut basis = Vec::with_capacity(d);
        for i in 0..d {
            let mut e = vec![0u32; d];
            e[i] = 1;
            basis.push(e);
        }
        basis
    }

    /// Creates a matrix from raw 2D row-major data, reducing entries modulo `p`.
    pub fn from_data(data: Vec<Vec<u32>>, p: u32) -> Self {
        let rows = data.len();
        let cols = if rows > 0 { data[0].len() } else { 0 };
        let mut normalized = vec![vec![0; cols]; rows];
        for r in 0..rows {
            for c in 0..cols {
                normalized[r][c] = data[r][c] % p;
            }
        }
        Self {
            data,
            p,
            rows,
            cols,
        }
    }

    /// Matrix addition modulo `p`.
    pub fn add(&self, rhs: &ModularMatrix) -> ModularMatrix {
        assert_eq!(self.rows, rhs.rows, "Matrix dimension mismatch in add");
        assert_eq!(self.cols, rhs.cols, "Matrix dimension mismatch in add");
        assert_eq!(self.p, rhs.p, "Field characteristic mismatch in add");

        let p = self.p;
        let mut data = vec![vec![0; self.cols]; self.rows];
        for r in 0..self.rows {
            for c in 0..self.cols {
                data[r][c] = (self.data[r][c] + rhs.data[r][c]) % p;
            }
        }
        ModularMatrix {
            data,
            p,
            rows: self.rows,
            cols: self.cols,
        }
    }

    /// Matrix subtraction modulo `p`.
    pub fn sub(&self, rhs: &ModularMatrix) -> ModularMatrix {
        assert_eq!(self.rows, rhs.rows, "Matrix dimension mismatch in sub");
        assert_eq!(self.cols, rhs.cols, "Matrix dimension mismatch in sub");
        assert_eq!(self.p, rhs.p, "Field characteristic mismatch in sub");

        let p = self.p;
        let mut data = vec![vec![0; self.cols]; self.rows];
        for r in 0..self.rows {
            for c in 0..self.cols {
                let diff = (self.data[r][c] + p - (rhs.data[r][c] % p)) % p;
                data[r][c] = diff;
            }
        }
        ModularMatrix {
            data,
            p,
            rows: self.rows,
            cols: self.cols,
        }
    }

    /// Matrix multiplication modulo `p`.
    pub fn mul(&self, rhs: &ModularMatrix) -> ModularMatrix {
        assert_eq!(self.cols, rhs.rows, "Matrix dimension mismatch in mul");
        assert_eq!(self.p, rhs.p, "Field characteristic mismatch in mul");

        let p = self.p as u64;
        let mut data = vec![vec![0; rhs.cols]; self.rows];
        for r in 0..self.rows {
            for c in 0..rhs.cols {
                let mut sum = 0u64;
                for k in 0..self.cols {
                    sum = (sum + (self.data[r][k] as u64 * rhs.data[k][c] as u64)) % p;
                }
                data[r][c] = sum as u32;
            }
        }
        ModularMatrix {
            data,
            p: self.p,
            rows: self.rows,
            cols: rhs.cols,
        }
    }

    /// Multiplies the matrix by a scalar `scalar` in `F_p`.
    pub fn scalar_mul(&self, scalar: u32) -> ModularMatrix {
        let p = self.p as u64;
        let s = (scalar % self.p) as u64;
        let mut data = vec![vec![0; self.cols]; self.rows];
        for r in 0..self.rows {
            for c in 0..self.cols {
                data[r][c] = ((self.data[r][c] as u64 * s) % p) as u32;
            }
        }
        ModularMatrix {
            data,
            p: self.p,
            rows: self.rows,
            cols: self.cols,
        }
    }

    /// Multiplies the matrix by a column vector `v` in `F_p^{cols}`.
    pub fn mul_vec(&self, v: &[u32]) -> Vec<u32> {
        assert_eq!(self.cols, v.len(), "Vector dimension mismatch in mul_vec");
        let p = self.p as u64;
        let mut result = vec![0u32; self.rows];
        for r in 0..self.rows {
            let mut sum = 0u64;
            for c in 0..self.cols {
                sum = (sum + (self.data[r][c] as u64 * (v[c] % self.p) as u64)) % p;
            }
            result[r] = sum as u32;
        }
        result
    }

    /// Computes the Reduced Row Echelon Form (RREF) in place using Gaussian elimination.
    ///
    /// Returns the indices of the pivot columns and the rank of the matrix.
    pub fn rref(&mut self) -> (Vec<usize>, usize) {
        let p = self.p;
        let mut pivot_row = 0;
        let mut pivot_cols = Vec::new();

        for col in 0..self.cols {
            if pivot_row >= self.rows {
                break;
            }

            // Find pivot in current column
            let mut cand_row = None;
            for r in pivot_row..self.rows {
                if self.data[r][col] % p != 0 {
                    cand_row = Some(r);
                    break;
                }
            }

            let Some(r) = cand_row else {
                continue;
            };

            // Swap rows
            self.data.swap(pivot_row, r);

            // Normalize pivot row
            let pivot_val = self.data[pivot_row][col] % p;
            let inv = mod_inverse(pivot_val, p).expect("pivot must be invertible");
            for c in 0..self.cols {
                self.data[pivot_row][c] =
                    ((self.data[pivot_row][c] as u64 * inv as u64) % p as u64) as u32;
            }

            // Eliminate all other rows
            for r in 0..self.rows {
                if r != pivot_row {
                    let factor = self.data[r][col] % p;
                    if factor != 0 {
                        for c in 0..self.cols {
                            let sub = (factor as u64 * self.data[pivot_row][c] as u64) % p as u64;
                            let cur = self.data[r][c] as u64;
                            self.data[r][c] = ((cur + p as u64 - sub) % p as u64) as u32;
                        }
                    }
                }
            }

            pivot_cols.push(col);
            pivot_row += 1;
        }

        (pivot_cols, pivot_row)
    }

    /// Computes a basis for the nullspace (kernel) of the matrix.
    pub fn nullspace(&self) -> Vec<Vec<u32>> {
        if self.cols == 0 {
            return Vec::new();
        }
        let p = self.p;
        let mut rref_mat = self.clone();
        let (pivot_cols, _) = rref_mat.rref();

        let mut is_pivot = vec![false; self.cols];
        for &c in &pivot_cols {
            is_pivot[c] = true;
        }

        let free_cols: Vec<usize> = (0..self.cols).filter(|&c| !is_pivot[c]).collect();
        let mut null_basis = Vec::with_capacity(free_cols.len());

        for &f in &free_cols {
            let mut v = vec![0u32; self.cols];
            v[f] = 1;

            for (row_idx, &p_col) in pivot_cols.iter().enumerate() {
                let coeff = rref_mat.data[row_idx][f] % p;
                v[p_col] = (p - coeff) % p;
            }

            null_basis.push(v);
        }

        null_basis
    }
}

/// A subspace of `F_p^d`, represented by a basis in row-echelon form.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModularSubspace {
    pub basis: Vec<Vec<u32>>,
    pub d: usize,
    pub p: u32,
}

impl ModularSubspace {
    /// Creates an empty subspace `{0}`.
    pub fn zero(d: usize, p: u32) -> Self {
        Self {
            basis: Vec::new(),
            d,
            p,
        }
    }

    /// Creates the full subspace `F_p^d`.
    pub fn full(d: usize, p: u32) -> Self {
        Self {
            basis: ModularMatrix::standard_basis(d),
            d,
            p,
        }
    }

    /// Adds a vector to the subspace and returns `true` if the dimension increased.
    pub fn add_vector(&mut self, mut v: Vec<u32>) -> bool {
        let p = self.p;
        for b in &self.basis {
            let lead_idx = b.iter().position(|&x| x != 0).unwrap();
            let factor = v[lead_idx];
            if factor != 0 {
                for i in lead_idx..self.d {
                    let sub = (factor as u64 * b[i] as u64) % p as u64;
                    v[i] = ((v[i] as u64 + p as u64 - sub) % p as u64) as u32;
                }
            }
        }

        if let Some(lead_idx) = v.iter().position(|&x| x != 0) {
            let inv = mod_inverse(v[lead_idx], p).unwrap();
            for i in lead_idx..self.d {
                v[i] = ((v[i] as u64 * inv as u64) % p as u64) as u32;
            }

            let insert_idx = self
                .basis
                .partition_point(|b| b.iter().position(|&x| x != 0).unwrap() < lead_idx);

            for b in &mut self.basis[..insert_idx] {
                let factor = b[lead_idx];
                if factor != 0 {
                    for i in lead_idx..self.d {
                        let sub = (factor as u64 * v[i] as u64) % p as u64;
                        b[i] = ((b[i] as u64 + p as u64 - sub) % p as u64) as u32;
                    }
                }
            }

            self.basis.insert(insert_idx, v);
            true
        } else {
            false
        }
    }

    /// Dimension of the subspace.
    pub fn dim(&self) -> usize {
        self.basis.len()
    }
}

/// Recursively computes a composition series of the `F_p[G]`-module `F_p^d`.
/// Returns a chain of submodules `0 = V_0 < V_1 < ... < V_k = V`.
pub fn composition_series(matrices: &[ModularMatrix], p: u32, d: usize) -> Vec<ModularSubspace> {
    if d == 0 {
        return vec![ModularSubspace::zero(d, p)];
    }

    let mut chain = vec![ModularSubspace::zero(d, p), ModularSubspace::full(d, p)];

    let mut refined = true;
    while refined {
        refined = false;
        let mut new_chain = Vec::new();

        for i in 0..chain.len() - 1 {
            let lower = &chain[i];
            let upper = &chain[i + 1];
            new_chain.push(lower.clone());

            if upper.dim() - lower.dim() > 1 {
                let mut found_intermediate = false;

                for b in &upper.basis {
                    let mut test = lower.clone();
                    if test.add_vector(b.clone()) {
                        let mut head = 0;
                        while head < test.dim() {
                            let cv = test.basis[head].clone();
                            for m in matrices {
                                test.add_vector(m.mul_vec(&cv));
                            }
                            head += 1;
                        }

                        if test.dim() < upper.dim() {
                            new_chain.push(test);
                            found_intermediate = true;
                            refined = true;
                            break;
                        }
                    }
                }

                if !found_intermediate {
                    'nullspace_search: for m in matrices {
                        for lambda in 0..p {
                            let lambda_id = ModularMatrix::identity(d, p).scalar_mul(lambda);
                            let shifted = m.sub(&lambda_id);
                            for null_v in shifted.nullspace() {
                                let mut test = lower.clone();
                                if test.add_vector(null_v) {
                                    let mut head = 0;
                                    while head < test.dim() {
                                        let cv = test.basis[head].clone();
                                        for m2 in matrices {
                                            test.add_vector(m2.mul_vec(&cv));
                                        }
                                        head += 1;
                                    }

                                    if upper.dim() == d && test.dim() < d {
                                        new_chain.push(test);
                                        found_intermediate = true;
                                        refined = true;
                                        break 'nullspace_search;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        new_chain.push(chain.last().unwrap().clone());
        chain = new_chain;
    }

    chain
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matrix_rref_nullspace() {
        let mut m = ModularMatrix::from_data(vec![vec![1, 2], vec![2, 1]], 3);
        let (pivots, rank) = m.rref();
        assert_eq!(pivots, vec![0]);
        assert_eq!(rank, 1);

        let ns = ModularMatrix::from_data(vec![vec![1, 2], vec![2, 1]], 3).nullspace();
        assert_eq!(ns.len(), 1);
        assert_eq!(ns[0], vec![1, 1]);
    }

    #[test]
    fn test_subspace_add() {
        let mut sub = ModularSubspace::zero(3, 2);
        assert!(sub.add_vector(vec![0, 1, 1]));
        assert!(sub.add_vector(vec![1, 1, 0]));
        assert!(!sub.add_vector(vec![1, 0, 1]));
        assert_eq!(sub.dim(), 2);
        assert_eq!(sub.basis, vec![vec![1, 0, 1], vec![0, 1, 1]]);
    }
}
