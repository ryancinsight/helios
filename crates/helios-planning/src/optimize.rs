//! Projected-gradient beam-weight optimization.

use helios_core::HeliosError;
use helios_math::{NumericElement, Scalar};

/// A dense linear dose-influence matrix `A` (rows = voxels, columns = beamlets):
/// `dose = A · x`. Row-major.
#[derive(Debug, Clone, PartialEq)]
pub struct DoseInfluence<T: Scalar> {
    voxels: usize,
    beamlets: usize,
    data: Vec<T>,
}

impl<T: Scalar> DoseInfluence<T> {
    /// Construct from a row-major `voxels × beamlets` matrix.
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if `data.len() != voxels·beamlets`.
    pub fn from_rows(voxels: usize, beamlets: usize, data: Vec<T>) -> Result<Self, HeliosError> {
        if data.len() != voxels * beamlets {
            return Err(HeliosError::InvalidDomainValue {
                field: "DoseInfluence::len",
                value: data.len() as f64,
                reason: "matrix length must equal voxels·beamlets",
            });
        }
        Ok(Self {
            voxels,
            beamlets,
            data,
        })
    }

    /// `(voxels, beamlets)`.
    #[must_use]
    pub fn dims(&self) -> (usize, usize) {
        (self.voxels, self.beamlets)
    }

    /// Zero-copy view of the row-major (`voxels × beamlets`) matrix entries.
    #[must_use]
    pub fn rows(&self) -> &[T] {
        &self.data
    }

    /// Dose `A · x` from beamlet weights `x` (length `beamlets`).
    #[must_use]
    pub fn apply(&self, x: &[T]) -> Vec<T> {
        let zero = <T as NumericElement>::ZERO;
        (0..self.voxels)
            .map(|i| {
                let row = &self.data[i * self.beamlets..(i + 1) * self.beamlets];
                row.iter().zip(x).fold(zero, |acc, (&a, &xj)| acc + a * xj)
            })
            .collect()
    }

    /// `Aᵀ · r` from a voxel-space residual `r` (length `voxels`).
    #[must_use]
    pub fn transpose_apply(&self, r: &[T]) -> Vec<T> {
        let zero = <T as NumericElement>::ZERO;
        let mut out = vec![zero; self.beamlets];
        for (i, &ri) in r.iter().enumerate() {
            let row = &self.data[i * self.beamlets..(i + 1) * self.beamlets];
            for (o, &a) in out.iter_mut().zip(row) {
                *o += a * ri;
            }
        }
        out
    }
}

/// Quadratic objective `½‖A x − d‖²`.
#[must_use]
pub fn objective_value<T: Scalar>(influence: &DoseInfluence<T>, x: &[T], prescription: &[T]) -> T {
    let zero = <T as NumericElement>::ZERO;
    let dose = influence.apply(x);
    let sum_sq = dose.iter().zip(prescription).fold(zero, |acc, (&di, &pi)| {
        let e = di - pi;
        acc + e * e
    });
    sum_sq * T::from_f64(0.5)
}

/// Minimize `½‖A x − d‖²` over `x ≥ 0` by projected gradient descent.
///
/// Iterates `x ← max(0, x − step·Aᵀ(A x − d))` from `x = 0`. For convergence the
/// `step` must satisfy `step < 2/‖AᵀA‖`. Returns the optimized non-negative
/// beamlet weights (length `beamlets`).
#[must_use]
pub fn optimize_beam_weights<T: Scalar>(
    influence: &DoseInfluence<T>,
    prescription: &[T],
    iterations: usize,
    step: T,
) -> Vec<T> {
    let zero = <T as NumericElement>::ZERO;
    let (_voxels, beamlets) = influence.dims();
    let mut x = vec![zero; beamlets];
    for _ in 0..iterations {
        let dose = influence.apply(&x);
        let residual: Vec<T> = dose
            .iter()
            .zip(prescription)
            .map(|(&di, &pi)| di - pi)
            .collect();
        let grad = influence.transpose_apply(&residual);
        for (xj, &gj) in x.iter_mut().zip(&grad) {
            *xj = (*xj - step * gj).max_scalar(zero);
        }
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    /// Identity influence: `n` voxels, `n` beamlets, `A = I`.
    fn identity(n: usize) -> DoseInfluence<f64> {
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            data[i * n + i] = 1.0;
        }
        DoseInfluence::from_rows(n, n, data).unwrap()
    }

    #[test]
    fn rejects_wrong_matrix_length() {
        assert!(DoseInfluence::from_rows(2, 3, vec![0.0; 5]).is_err());
        assert!(DoseInfluence::from_rows(2, 3, vec![0.0; 6]).is_ok());
    }

    #[test]
    fn apply_and_transpose_are_consistent() {
        // A = [[1,2],[0,1],[3,0]] (3 voxels, 2 beamlets).
        let a = DoseInfluence::from_rows(3, 2, vec![1.0, 2.0, 0.0, 1.0, 3.0, 0.0]).unwrap();
        assert_eq!(a.apply(&[1.0, 1.0]), vec![3.0, 1.0, 3.0]);
        // Aᵀ·[1,1,1] = column sums = [1+0+3, 2+1+0] = [4, 3].
        assert_eq!(a.transpose_apply(&[1.0, 1.0, 1.0]), vec![4.0, 3.0]);
    }

    #[test]
    fn identity_problem_converges_to_prescription() {
        // min ½‖x − d‖² s.t. x≥0, d≥0 → x = d.
        let a = identity(3);
        let d = [1.0, 2.0, 3.0];
        let x = optimize_beam_weights(&a, &d, 500, 0.5);
        for (xi, di) in x.iter().zip(&d) {
            assert_relative_eq!(xi, di, epsilon = 1e-6);
        }
    }

    #[test]
    fn negative_target_is_clamped_to_zero() {
        // A=I, prescription −1 for a voxel → unconstrained min is −1, projected to 0.
        let a = identity(2);
        let x = optimize_beam_weights(&a, &[-1.0, 2.0], 500, 0.5);
        assert_relative_eq!(x[0], 0.0, epsilon = 1e-9);
        assert_relative_eq!(x[1], 2.0, epsilon = 1e-6);
    }

    #[test]
    fn objective_decreases_monotonically() {
        // Over-determined: A = [[1,1],[1,-1],[2,0]], d = [2,0,2].
        let a = DoseInfluence::from_rows(3, 2, vec![1.0, 1.0, 1.0, -1.0, 2.0, 0.0]).unwrap();
        let d = [2.0, 0.0, 2.0];
        let mut x = vec![0.0; 2];
        let mut prev = objective_value(&a, &x, &d);
        for _ in 0..50 {
            x = optimize_beam_weights(&a, &d, 1, 0.1);
            let cur = objective_value(&a, &x, &d);
            assert!(cur <= prev + 1e-12, "objective rose: {prev} → {cur}");
            prev = cur;
        }
    }

    #[test]
    fn least_squares_solution_for_well_conditioned_problem() {
        // A = diag(2, 4); min ½‖A x − d‖² → x = d/diag (all positive).
        let a = DoseInfluence::from_rows(2, 2, vec![2.0, 0.0, 0.0, 4.0]).unwrap();
        let d = [6.0, 8.0];
        let x = optimize_beam_weights(&a, &d, 2000, 0.05);
        assert_relative_eq!(x[0], 3.0, epsilon = 1e-4); // 6/2
        assert_relative_eq!(x[1], 2.0, epsilon = 1e-4); // 8/4
    }

    #[test]
    fn optimizer_is_generic_over_scalar_f32() {
        let mut data = vec![0.0_f32; 4];
        data[0] = 1.0;
        data[3] = 1.0;
        let a = DoseInfluence::from_rows(2, 2, data).unwrap();
        let x = optimize_beam_weights(&a, &[1.5_f32, 0.5], 300, 0.5);
        assert_relative_eq!(x[0], 1.5_f32, epsilon = 1e-3);
        assert_relative_eq!(x[1], 0.5_f32, epsilon = 1e-3);
    }
}
