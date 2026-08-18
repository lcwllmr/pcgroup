//! Construction of irreducible representations of supersolvable groups via the Baum-Clausen algorithm.
//!
//! For finite supersolvable groups given by a power-commutator presentation, all irreducible
//! representations adapted to a chief series are monomial over roots of unity. The algorithm
//! operates bottom-up along the chief series in `O(|G| log |G|)` basic operations in `Z/eZ`.

#![allow(clippy::needless_range_loop)]

use std::collections::{HashMap, HashSet};
use std::fmt;

use crate::util::{MonomialMatrix, extract_roots, lcm};
use crate::{Element, GeneratingSequence, Presentation, Word, chief_series, is_supersolvable};

/// Errors arising during the Baum-Clausen representation induction algorithm.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaumClausenError {
    /// The given group presentation is not supersolvable.
    NotSupersolvable,
    /// Missing matrix for a generator in a representation.
    MissingMatrix {
        generator: usize,
        representation_id: usize,
    },
    /// A referenced representation ID was not found.
    MissingRepresentation { id: usize },
    /// Missing tau (conjugation tracking) map entry.
    MissingTauEntry {
        generator: usize,
        representation_id: usize,
    },
    /// Missing X (intertwining matrix) entry.
    MissingXEntry {
        generator: usize,
        representation_id: usize,
    },
    /// Root extraction failed when computing `p`-th roots of unity.
    RootExtractionFailed {
        shift: u32,
        order: u32,
        exponent: u32,
    },
    /// Induced representation has no children orbit recorded.
    MissingChildren { representation_id: usize },
    /// Invalid transition between stages in the chief series.
    InvalidStageTransition,
}

impl fmt::Display for BaumClausenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotSupersolvable => write!(f, "Group is not supersolvable."),
            Self::MissingMatrix {
                generator,
                representation_id,
            } => write!(
                f,
                "Missing matrix for generator g{generator} in representation {representation_id}."
            ),
            Self::MissingRepresentation { id } => {
                write!(f, "Representation with ID {id} not found.")
            }
            Self::MissingTauEntry {
                generator,
                representation_id,
            } => write!(
                f,
                "Missing tau entry for generator g{generator} and representation {representation_id}."
            ),
            Self::MissingXEntry {
                generator,
                representation_id,
            } => write!(
                f,
                "Missing X entry for generator g{generator} and representation {representation_id}."
            ),
            Self::RootExtractionFailed {
                shift,
                order,
                exponent,
            } => write!(
                f,
                "Failed to extract {order}-th root for shift {shift} modulo {exponent}."
            ),
            Self::MissingChildren { representation_id } => write!(
                f,
                "Induced representation {representation_id} has no children recorded."
            ),
            Self::InvalidStageTransition => write!(f, "Invalid stage transition encountered."),
        }
    }
}

impl std::error::Error for BaumClausenError {}

/// An irreducible monomial representation constructed at a stage of the Baum-Clausen algorithm.
///
/// # Examples
/// ```
/// use pcgroup::zoo::dihedral;
/// use pcgroup::{Element, irreducible_representations};
///
/// let s3 = dihedral(3);
/// let irreps = irreducible_representations(&s3).unwrap();
///
/// // S_3 has two 1D irreps and one 2D irrep
/// let rep2d = &irreps[2];
/// assert_eq!(rep2d.dim, 2);
///
/// // Evaluating reflection g0: has order 2
/// let g0 = Element::from_generator(0, 1, &s3);
/// let m_g0 = rep2d.evaluate_element(&g0);
/// assert_eq!(m_g0.pow(2), rep2d.evaluate_element(&Element::identity(&s3)));
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Representation {
    /// Unique sequentially assigned representation identifier.
    pub id: usize,
    /// Dimension of the representation vector space.
    pub dim: usize,
    /// Root-of-unity exponent `e`.
    pub exponent: u32,
    /// Monomial matrix for each generator `g0, ..., g{n-1}`.
    pub matrices: Vec<MonomialMatrix>,
    /// Whether this representation was constructed via induction (`true`) or extension (`false`).
    pub is_induced: bool,
    /// IDs of child orbit representations from the previous stage that induced this representation.
    pub children: Option<Vec<usize>>,
    /// For Case 1 (extension), IDs of sibling extensions sharing the same base representation.
    pub cluster: Option<Vec<usize>>,
}

impl Representation {
    /// Evaluates a group element given as a [`Word`] to its [`MonomialMatrix`] in this representation.
    pub fn evaluate(&self, word: &Word) -> MonomialMatrix {
        let mut curr = MonomialMatrix::identity(self.dim, self.exponent);
        let mut next = MonomialMatrix::identity(self.dim, self.exponent);
        for term in word.iter() {
            let matrix = &self.matrices[term.gen_index];
            let factor = matrix.pow(term.exponent as i64);
            curr.compose_into(&factor, &mut next);
            std::mem::swap(&mut curr, &mut next);
        }
        curr
    }

    /// Evaluates a group [`Element`] to its [`MonomialMatrix`] in this representation.
    #[inline]
    pub fn evaluate_element(&self, element: &Element<'_>) -> MonomialMatrix {
        self.evaluate(element.word())
    }
}

/// Represents a stage in the Baum-Clausen induction algorithm along a chief series.
///
/// # Examples
/// ```
/// use pcgroup::zoo::dihedral;
/// use pcgroup::{BaumClausenStage, chief_series};
///
/// let s3 = dihedral(3);
/// let chief = chief_series(&s3);
///
/// // Start at level 0 with the trivial group {1}
/// let mut stage = BaumClausenStage::trivial(&s3, 6);
/// assert_eq!(stage.level, 0);
/// assert_eq!(stage.representations.len(), 1);
///
/// // Advance bottom-up along the chief series: {1} < C_3 < S_3
/// for next_subgroup in chief.iter().rev().skip(1) {
///     stage = stage.next(next_subgroup).unwrap();
/// }
/// assert_eq!(stage.representations.len(), 3);
/// ```
#[derive(Clone, Debug)]
pub struct BaumClausenStage<'a> {
    /// Current induction level (step index in the chief series, 0 for `{1}`).
    pub level: usize,
    /// Shared root-of-unity exponent `e` for all monomial matrices.
    pub exponent: u32,
    /// Irreducible representations computed at this stage.
    pub representations: Vec<Representation>,
    /// Intertwining matrices `X(g_j, F)` for higher generators `g_j` acting on representations `F`.
    pub x_dict: HashMap<usize, HashMap<usize, MonomialMatrix>>,
    /// Conjugation tracking map `tau(g_j, F)` mapping a representation to its conjugate representation ID.
    pub tau_dict: HashMap<usize, HashMap<usize, usize>>,
    /// Unprocessed generators belonging to higher levels in the chief series.
    pub higher_generators: Vec<usize>,
    /// Associated presentation.
    pub pres: &'a Presentation,
}

impl<'a> BaumClausenStage<'a> {
    /// Initializes the trivial stage (level 0) with a single 1-dimensional trivial representation.
    pub fn trivial(pres: &'a Presentation, exponent: u32) -> Self {
        let n_gens = pres.num_gens();
        let matrices = vec![MonomialMatrix::identity(1, exponent); n_gens];

        let trivial_rep = Representation {
            id: 0,
            dim: 1,
            exponent,
            matrices,
            is_induced: false,
            children: None,
            cluster: None,
        };

        let generators: Vec<usize> = (0..n_gens).collect();
        let mut x_dict = HashMap::new();
        let mut tau_dict = HashMap::new();

        for &g in &generators {
            let mut x_map = HashMap::new();
            x_map.insert(0, MonomialMatrix::identity(1, exponent));
            x_dict.insert(g, x_map);

            let mut tau_map = HashMap::new();
            tau_map.insert(0, 0);
            tau_dict.insert(g, tau_map);
        }

        Self {
            level: 0,
            exponent,
            representations: vec![trivial_rep],
            x_dict,
            tau_dict,
            higher_generators: generators,
            pres,
        }
    }

    /// Advances the algorithm from the current stage to the next level in the chief series.
    pub fn next(&self, next_subgroup: &GeneratingSequence<'a>) -> Result<Self, BaumClausenError> {
        let mut found_gen = None;
        for elem in next_subgroup.elements() {
            let lead = elem.leading_generator().unwrap();
            if self.higher_generators.contains(&lead) {
                found_gen = Some(lead);
                break;
            }
        }

        let Some(g_i) = found_gen else {
            return Err(BaumClausenError::InvalidStageTransition);
        };

        let p = self.pres.relative_order(g_i);
        self.advance_step(g_i, p)
    }

    /// Advances along generator `g_i` of relative order `p`.
    pub(crate) fn advance_step(&self, g_i: usize, p: u32) -> Result<Self, BaumClausenError> {
        let level = self.level + 1;
        let e = self.exponent;
        let n_gens = self.pres.num_gens();
        let higher_gens: Vec<usize> = self
            .higher_generators
            .iter()
            .copied()
            .filter(|&g| g != g_i)
            .collect();

        let id_to_rep: HashMap<usize, &Representation> = self
            .representations
            .iter()
            .map(|rep| (rep.id, rep))
            .collect();

        let mut visited = HashSet::new();
        let mut new_representations = Vec::new();
        let mut next_id = 0;

        let mut case1_clusters: Vec<(usize, Vec<usize>)> = Vec::new();
        let mut case2_orbits: Vec<(Representation, Vec<usize>)> = Vec::new();

        let power_word_gi = self.pres.power(g_i);

        // ==========================================
        // PHASE 1: Induce / Extend Representations
        // ==========================================
        for f in &self.representations {
            if visited.contains(&f.id) {
                continue;
            }

            let tau_gi = *self.tau_dict.get(&g_i).and_then(|m| m.get(&f.id)).ok_or(
                BaumClausenError::MissingTauEntry {
                    generator: g_i,
                    representation_id: f.id,
                },
            )?;

            if tau_gi == f.id {
                // Case 1: Invariant representation under g_i (Extension)
                visited.insert(f.id);

                let x_if = self.x_dict.get(&g_i).and_then(|m| m.get(&f.id)).ok_or(
                    BaumClausenError::MissingXEntry {
                        generator: g_i,
                        representation_id: f.id,
                    },
                )?;

                let f_gp = f.evaluate(power_word_gi);
                let x_if_p = x_if.pow(p as i64);

                let target_val = f_gp.vals.first().copied().unwrap_or(0);
                let current_val = x_if_p.vals.first().copied().unwrap_or(0);
                let shift = (target_val + e - (current_val % e)) % e;
                let roots = extract_roots(shift, p, e);

                if roots.len() < p as usize {
                    return Err(BaumClausenError::RootExtractionFailed {
                        shift,
                        order: p,
                        exponent: e,
                    });
                }

                let cluster_ids: Vec<usize> = (next_id..(next_id + (p as usize))).collect();
                let mut extensions = Vec::with_capacity(p as usize);

                for c_k in roots.into_iter().take(p as usize) {
                    let mut matrices = f.matrices.clone();
                    matrices[g_i] = x_if.scalar_mul(c_k);

                    let d_k = Representation {
                        id: next_id,
                        dim: f.dim,
                        exponent: e,
                        matrices,
                        is_induced: false,
                        children: Some(vec![f.id]),
                        cluster: Some(cluster_ids.clone()),
                    };
                    next_id += 1;
                    extensions.push(d_k);
                }

                case1_clusters.push((f.id, cluster_ids));
                new_representations.extend(extensions);
            } else {
                // Case 2: Orbit of size p under g_i (Induction)
                let mut orbit: Vec<&Representation> = vec![f];
                visited.insert(f.id);
                let mut cursor = f;

                for _ in 0..(p - 1) {
                    let next_orbit_id = *self
                        .tau_dict
                        .get(&g_i)
                        .and_then(|m| m.get(&cursor.id))
                        .ok_or(BaumClausenError::MissingTauEntry {
                            generator: g_i,
                            representation_id: cursor.id,
                        })?;
                    cursor = id_to_rep
                        .get(&next_orbit_id)
                        .copied()
                        .ok_or(BaumClausenError::MissingRepresentation { id: next_orbit_id })?;
                    orbit.push(cursor);
                    visited.insert(cursor.id);
                }

                let dim = f.dim;
                let full_dim = (p as usize) * dim;
                let mut matrices = vec![MonomialMatrix::identity(full_dim, e); n_gens];

                // Direct sums for lower generators
                for g in 0..n_gens {
                    if g != g_i {
                        let mut acc = MonomialMatrix::new(Vec::new(), Vec::new(), e);
                        for rep in &orbit {
                            acc = acc.direct_sum(&rep.matrices[g]);
                        }
                        matrices[g] = acc;
                    }
                }

                // Block intertwiners X_k
                let mut x_blocks = vec![MonomialMatrix::identity(dim, e)];
                for k in 1..(p as usize) {
                    let x_prev = &x_blocks[k - 1];
                    let x_k = self
                        .x_dict
                        .get(&g_i)
                        .and_then(|m| m.get(&orbit[k - 1].id))
                        .ok_or(BaumClausenError::MissingXEntry {
                            generator: g_i,
                            representation_id: orbit[k - 1].id,
                        })?
                        .compose(x_prev);
                    x_blocks.push(x_k);
                }

                let f0_gp = orbit[0].evaluate(power_word_gi);
                let b_last = f0_gp.compose(&x_blocks[(p as usize) - 1].invert());

                let mut gi_perm = vec![0; full_dim];
                let mut gi_vals = vec![0; full_dim];

                for k in 0..(p as usize) {
                    let (block_mat, target_block) = if k < (p as usize) - 1 {
                        (
                            self.x_dict
                                .get(&g_i)
                                .and_then(|m| m.get(&orbit[k].id))
                                .unwrap(),
                            k + 1,
                        )
                    } else {
                        (&b_last, 0)
                    };

                    for r in 0..dim {
                        let global_row = k * dim + r;
                        let global_col = target_block * dim + block_mat.perm[r];
                        gi_perm[global_row] = global_col;
                        gi_vals[global_row] = block_mat.vals[r];
                    }
                }

                matrices[g_i] = MonomialMatrix::new(gi_perm, gi_vals, e);

                let orbit_ids: Vec<usize> = orbit.iter().map(|rep| rep.id).collect();
                let d = Representation {
                    id: next_id,
                    dim: full_dim,
                    exponent: e,
                    matrices,
                    is_induced: true,
                    children: Some(orbit_ids.clone()),
                    cluster: None,
                };
                next_id += 1;
                case2_orbits.push((d.clone(), orbit_ids));
                new_representations.push(d);
            }
        }

        // ==========================================
        // PHASE 2: Higher Intertwiners & Conjugations
        // ==========================================
        let mut x_dict: HashMap<usize, HashMap<usize, MonomialMatrix>> = higher_gens
            .iter()
            .copied()
            .map(|g_j| (g_j, HashMap::new()))
            .collect();
        let mut tau_dict: HashMap<usize, HashMap<usize, usize>> = higher_gens
            .iter()
            .copied()
            .map(|g_j| (g_j, HashMap::new()))
            .collect();

        let new_id_to_rep: HashMap<usize, &Representation> = new_representations
            .iter()
            .map(|rep| (rep.id, rep))
            .collect();

        let rep_to_induced: HashMap<usize, &Representation> = case2_orbits
            .iter()
            .flat_map(|(d, orbit_ids)| orbit_ids.iter().map(move |&rep_id| (rep_id, d)))
            .collect();

        let base_to_cluster: HashMap<usize, &[usize]> = case1_clusters
            .iter()
            .map(|(base_id, cluster_ids)| (*base_id, cluster_ids.as_slice()))
            .collect();

        for &g_j in &higher_gens {
            let g_j_elem = Element::from_generator(g_j, 1, self.pres);
            let g_i_elem = Element::from_generator(g_i, 1, self.pres);
            let conj_elem = &g_j_elem.inverse() * &(&g_i_elem * &g_j_elem);
            let a_j = conj_elem
                .word()
                .iter()
                .find(|t| t.gen_index == g_i)
                .map_or(1, |t| t.exponent);

            // Phase 2 Case 1: Extensions
            for &(base_id, ref extension_ids) in &case1_clusters {
                let pi_j_base_id = *self
                    .tau_dict
                    .get(&g_j)
                    .and_then(|m| m.get(&base_id))
                    .ok_or(BaumClausenError::MissingTauEntry {
                        generator: g_j,
                        representation_id: base_id,
                    })?;

                let target_cluster_ids = base_to_cluster
                    .get(&pi_j_base_id)
                    .copied()
                    .ok_or(BaumClausenError::MissingRepresentation { id: pi_j_base_id })?;

                let x_jf = self.x_dict.get(&g_j).and_then(|m| m.get(&base_id)).ok_or(
                    BaumClausenError::MissingXEntry {
                        generator: g_j,
                        representation_id: base_id,
                    },
                )?;

                let d0 = new_id_to_rep[&extension_ids[0]];
                let d0_w = d0.evaluate_element(&conj_elem);
                let m_mat = x_jf.compose(&d0_w).compose(&x_jf.invert());

                let target_phase = m_mat.vals.first().copied().unwrap_or(0);

                let mut found_l = 0;
                for (l, &cand_id) in target_cluster_ids.iter().enumerate() {
                    let cand_rep = new_id_to_rep[&cand_id];
                    let cand_phase = cand_rep.matrices[g_i].vals.first().copied().unwrap_or(0);
                    if cand_phase == target_phase {
                        found_l = l;
                        break;
                    }
                }

                let p_exts = extension_ids.len();
                for (k, &d_k_id) in extension_ids.iter().enumerate() {
                    let target_idx = (found_l + k * (a_j as usize)) % p_exts;
                    let target_id = target_cluster_ids[target_idx];

                    tau_dict.get_mut(&g_j).unwrap().insert(d_k_id, target_id);
                    x_dict.get_mut(&g_j).unwrap().insert(d_k_id, x_jf.clone());
                }
            }

            // Phase 2 Case 2: Induced
            for (d, orbit_ids) in &case2_orbits {
                let f0_id = orbit_ids[0];
                let f0_dim = id_to_rep[&f0_id].dim;

                let pi_j_f0_id = *self.tau_dict.get(&g_j).and_then(|m| m.get(&f0_id)).ok_or(
                    BaumClausenError::MissingTauEntry {
                        generator: g_j,
                        representation_id: f0_id,
                    },
                )?;

                let target_d = rep_to_induced
                    .get(&pi_j_f0_id)
                    .copied()
                    .ok_or(BaumClausenError::MissingRepresentation { id: pi_j_f0_id })?;

                tau_dict.get_mut(&g_j).unwrap().insert(d.id, target_d.id);

                let target_orbit_ids = target_d.children.as_ref().unwrap();

                let mut sigma = Vec::with_capacity(p as usize);
                let mut x_k_list = Vec::with_capacity(p as usize);

                for &fk_id in orbit_ids {
                    let pi_j_fk_id = *self.tau_dict[&g_j].get(&fk_id).unwrap();
                    let sig_k = target_orbit_ids
                        .iter()
                        .position(|&rep_id| rep_id == pi_j_fk_id)
                        .unwrap();
                    sigma.push(sig_k);

                    let x_k = self.x_dict[&g_j].get(&fk_id).unwrap().clone();
                    x_k_list.push(x_k);
                }

                let mut inv_sigma = vec![0; p as usize];
                for (k, &s) in sigma.iter().enumerate() {
                    inv_sigma[s] = k;
                }

                let d_w = d.evaluate_element(&conj_elem);
                let target_d_gi = &target_d.matrices[g_i];

                let mut c_phases = vec![0u32; p as usize];
                let mut cur_k = 0;

                for _ in 0..(p - 1) as usize {
                    let r = 0;
                    let v = sigma[cur_k] * f0_dim + r;
                    let u = cur_k * f0_dim + x_k_list[cur_k].perm[r];
                    let v_prime = target_d_gi.perm[v];
                    let target_block = v_prime / f0_dim;
                    let k_prev = inv_sigma[target_block];
                    let r_prime = v_prime % f0_dim;

                    let lhs_term = (x_k_list[cur_k].vals[r] + d_w.vals[u]) % e;
                    let rhs_term = (target_d_gi.vals[v] + x_k_list[k_prev].vals[r_prime]) % e;
                    let diff = (lhs_term + e - rhs_term) % e;

                    c_phases[k_prev] = (c_phases[cur_k] + diff) % e;
                    cur_k = k_prev;
                }

                let mut full_perm = vec![0; d.dim];
                let mut full_vals = vec![0; d.dim];

                for k in 0..(p as usize) {
                    let sig_k = sigma[k];
                    let x_k = &x_k_list[k];
                    let c_scalar = c_phases[k];

                    for r in 0..f0_dim {
                        let row = sig_k * f0_dim + r;
                        let col = k * f0_dim + x_k.perm[r];
                        full_perm[row] = col;
                        full_vals[row] = (x_k.vals[r] + c_scalar) % e;
                    }
                }

                let y_mat = MonomialMatrix::new(full_perm, full_vals, e);
                x_dict.get_mut(&g_j).unwrap().insert(d.id, y_mat);
            }
        }

        Ok(Self {
            level,
            exponent: e,
            representations: new_representations,
            x_dict,
            tau_dict,
            higher_generators: higher_gens,
            pres: self.pres,
        })
    }
}

/// Computes the group exponent `e = lcm(ord(g))` for a polycyclic presentation.
pub fn group_exponent(pres: &Presentation) -> u32 {
    let n = pres.num_gens();
    if n == 0 {
        return 1;
    }

    let mut exp = 1u32;
    for i in 0..n {
        let elem = Element::from_generator(i, 1, pres);
        let mut cur = elem.clone();
        let mut ord = 1u32;
        while !cur.is_identity() {
            cur = &cur * &elem;
            ord += 1;
            if ord > 10000 {
                break;
            }
        }
        exp = lcm(exp, ord);
    }

    // Also include relative orders
    for &p in pres.relative_orders() {
        exp = lcm(exp, p);
    }

    exp
}

/// Computes a full set of pairwise inequivalent irreducible representations of a supersolvable group.
///
/// # Examples
/// ```
/// use pcgroup::zoo::{dihedral, quaternion};
/// use pcgroup::irreducible_representations;
///
/// // S_3 = D_6 of order 6 has 3 irreducible representations (two 1D, one 2D)
/// let s3 = dihedral(3);
/// let irreps_s3 = irreducible_representations(&s3).unwrap();
/// assert_eq!(irreps_s3.iter().map(|irrep| irrep.dim).collect::<Vec<_>>(), vec![1, 1, 2]);
///
/// // Q_8 of order 8 has 5 irreducible representations (four 1D, one 2D)
/// let q8 = quaternion(2);
/// let irreps_q8 = irreducible_representations(&q8).unwrap();
/// assert_eq!(irreps_q8.iter().map(|irrep| irrep.dim).collect::<Vec<_>>(), vec![1, 1, 1, 1, 2]);
/// ```
pub fn irreducible_representations(
    pres: &Presentation,
) -> Result<Vec<Representation>, BaumClausenError> {
    if !is_supersolvable(pres) {
        return Err(BaumClausenError::NotSupersolvable);
    }

    let chief = chief_series(pres);
    let exponent = group_exponent(pres);
    let mut stage = BaumClausenStage::trivial(pres, exponent);

    // Chief series is descending G = N_0 > N_1 > ... > N_m = {1}
    // We iterate bottom-up along the ascending chain {1} = G_0 < G_1 < ... < G_m = G
    for next_subgroup in chief.iter().rev().skip(1) {
        stage = stage.next(next_subgroup)?;
    }

    Ok(stage.representations)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zoo::{cyclic, dihedral, quaternion};

    #[test]
    fn test_error_display() {
        let err = BaumClausenError::MissingMatrix {
            generator: 1,
            representation_id: 3,
        };
        assert_eq!(
            err.to_string(),
            "Missing matrix for generator g1 in representation 3."
        );
    }

    #[test]
    fn test_trivial_stage() {
        let pres = cyclic(6);
        let stage = BaumClausenStage::trivial(&pres, 6);
        assert_eq!(stage.level, 0);
        assert_eq!(stage.representations.len(), 1);
        assert_eq!(stage.representations[0].dim, 1);
        assert_eq!(stage.higher_generators, vec![0, 1]);
    }

    fn verify_representation_homomorphism(pres: &Presentation, rep: &Representation) {
        // Check power relations: D(g_i)^{p_i} == D(power_word_i)
        for i in 0..pres.num_gens() {
            let p_i = pres.relative_order(i);
            let d_gi_p = rep.matrices[i].pow(p_i as i64);
            let d_power_word = rep.evaluate(pres.power(i));
            assert_eq!(
                d_gi_p, d_power_word,
                "Power relation failed for generator g{i} in representation ID {}",
                rep.id
            );
        }

        // Check commutator relations: [D(g_i), D(g_j)] == D(comm_word_ij)
        for i in 1..pres.num_gens() {
            for j in 0..i {
                let gi = &rep.matrices[i];
                let gj = &rep.matrices[j];
                // [gi, gj] = gi^-1 * gj^-1 * gi * gj
                let comm_matrix = gi.invert().compose(&gj.invert()).compose(gi).compose(gj);
                let d_comm_word = rep.evaluate(pres.commutator(i, j));
                assert_eq!(
                    comm_matrix, d_comm_word,
                    "Commutator relation failed for [g{i}, g{j}] in representation ID {}",
                    rep.id
                );
            }
        }
    }

    #[test]
    fn test_s3_irreps_homomorphisms() {
        let pres = dihedral(3); // S_3, order 6
        let irreps = irreducible_representations(&pres).expect("succeeds for S_3");
        assert_eq!(irreps.len(), 3);
        let dims: Vec<usize> = irreps.iter().map(|r| r.dim).collect();
        assert_eq!(dims, vec![1, 1, 2]);

        for rep in &irreps {
            verify_representation_homomorphism(&pres, rep);
        }
    }

    #[test]
    fn test_q8_irreps_homomorphisms() {
        let pres = quaternion(2); // Q_8, order 8
        let irreps = irreducible_representations(&pres).expect("succeeds for Q_8");
        assert_eq!(irreps.len(), 5);
        let dims: Vec<usize> = irreps.iter().map(|r| r.dim).collect();
        assert_eq!(dims, vec![1, 1, 1, 1, 2]);

        for rep in &irreps {
            verify_representation_homomorphism(&pres, rep);
        }
    }

    #[test]
    fn test_d8_irreps_homomorphisms() {
        let pres = dihedral(4); // D_8, order 8
        let irreps = irreducible_representations(&pres).expect("succeeds for D_8");
        assert_eq!(irreps.len(), 5);
        let dims: Vec<usize> = irreps.iter().map(|r| r.dim).collect();
        assert_eq!(dims, vec![1, 1, 1, 1, 2]);

        for rep in &irreps {
            verify_representation_homomorphism(&pres, rep);
        }
    }
}
