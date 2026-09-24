//! Projected-gradient beam-weight optimization.

use helios_core::HeliosError;
use helios_math::{NumericElement, Scalar};
use mnemosyne_arena::{AlignedVec, ScratchElement};

/// A dense linear dose-influence matrix `A` (rows = voxels, columns = beamlets):
/// `dose = A · x`. Row-major.
///
/// Entries are held in a 64-byte cache-line-aligned [`AlignedVec<T>`], so the
/// row-wise dot product in [`apply`](Self::apply) and the row-wise accumulation
/// in [`transpose_apply`](Self::transpose_apply) walk cache-line-aligned rows
/// without a realignment fixup.
///
/// # The `ScratchElement` bound
///
/// [`AlignedVec`] only admits [`ScratchElement`] elements, so `T` carries that
/// bound in addition to [`Scalar`]. This is not a narrowing of the accepted
/// set: [`Scalar`] is `eunomia::RealField`, whose supertrait `NumericElement` is
/// sealed and whose implementors are exactly `f32` and `f64` — both of which
/// implement `ScratchElement`.
#[derive(Debug, Clone, PartialEq)]
pub struct DoseInfluence<T: Scalar + ScratchElement> {
    voxels: usize,
    beamlets: usize,
    data: AlignedVec<T>,
}

impl<T: Scalar + ScratchElement> DoseInfluence<T> {
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
            // The caller hands the matrix over, so consume the `Vec` rather than
            // borrowing it. `AlignedVec` cannot adopt the allocation (the global
            // allocator gives no alignment guarantee), so this copies once.
            data: data.into(),
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
        let mut out = vec![<T as NumericElement>::ZERO; self.voxels];
        self.apply_into(x, &mut out);
        out
    }

    /// Writes `A · x` into `out` (length `voxels`), reusing the caller's buffer.
    ///
    /// The row-wise dot product is the innermost loop of the optimizer, so it is
    /// kept separate from [`apply`](Self::apply) to let the iteration reuse one
    /// dose buffer instead of allocating a fresh vector per step.
    fn apply_into(&self, x: &[T], out: &mut [T]) {
        let zero = <T as NumericElement>::ZERO;
        for (i, o) in out.iter_mut().enumerate() {
            let row = &self.data[i * self.beamlets..(i + 1) * self.beamlets];
            *o = row.iter().zip(x).fold(zero, |acc, (&a, &xj)| acc + a * xj);
        }
    }

    /// `Aᵀ · r` from a voxel-space residual `r` (length `voxels`).
    #[must_use]
    pub fn transpose_apply(&self, r: &[T]) -> Vec<T> {
        let mut out = vec![<T as NumericElement>::ZERO; self.beamlets];
        self.transpose_apply_into(r, &mut out);
        out
    }

    /// Writes `Aᵀ · r` into `out` (length `beamlets`), reusing the caller's
    /// buffer. `out` is overwritten, not accumulated into.
    ///
    /// Like [`apply_into`](Self::apply_into), this exists so the optimizer can
    /// keep one gradient buffer across all iterations.
    fn transpose_apply_into(&self, r: &[T], out: &mut [T]) {
        out.fill(<T as NumericElement>::ZERO);
        for (i, &ri) in r.iter().enumerate() {
            let row = &self.data[i * self.beamlets..(i + 1) * self.beamlets];
            for (o, &a) in out.iter_mut().zip(row) {
                *o += a * ri;
            }
        }
    }
}

/// Quadratic objective `½‖A x − d‖²`.
#[must_use]
pub fn objective_value<T: Scalar + ScratchElement>(
    influence: &DoseInfluence<T>,
    x: &[T],
    prescription: &[T],
) -> T {
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
///
/// The three voxel-/beamlet-space work vectors are allocated once and reused
/// across the iteration, so the loop body itself performs no allocation.
#[must_use]
pub fn optimize_beam_weights<T: Scalar + ScratchElement>(
    influence: &DoseInfluence<T>,
    prescription: &[T],
    iterations: usize,
    step: T,
) -> Vec<T> {
    let zero = <T as NumericElement>::ZERO;
    let (voxels, beamlets) = influence.dims();
    let mut x = vec![zero; beamlets];
    let mut dose = vec![zero; voxels];
    let mut residual = vec![zero; voxels];
    let mut grad = vec![zero; beamlets];
    for _ in 0..iterations {
        influence.apply_into(&x, &mut dose);
        for (r, (&di, &pi)) in residual.iter_mut().zip(dose.iter().zip(prescription)) {
            *r = di - pi;
        }
        influence.transpose_apply_into(&residual, &mut grad);
        for (xj, &gj) in x.iter_mut().zip(&grad) {
            *xj = (*xj - step * gj).max_scalar(zero);
        }
    }
    x
}

#[cfg(test)]
mod tests {
    #![expect(
        clippy::unwrap_used,
        reason = "ratchet HELIOS-UNWRAP-1: pre-existing debt"
    )]
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_math::ShippedScalar;

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
        let m = DoseInfluence::from_rows(2, 3, vec![0.0; 6]).unwrap();
        assert_eq!(m.dims(), (2, 3));
        assert_eq!(m.rows(), &[0.0; 6]);
    }

    #[test]
    fn apply_and_transpose_are_consistent() {
        // A = [[1,2],[0,1],[3,0]] (3 voxels, 2 beamlets).
        let a = DoseInfluence::from_rows(3, 2, vec![1.0, 2.0, 0.0, 1.0, 3.0, 0.0]).unwrap();
        assert_eq!(a.apply(&[1.0, 1.0]), vec![3.0, 1.0, 3.0]);
        // Aᵀ·[1,1,1] = column sums = [1+0+3, 2+1+0] = [4, 3].
        assert_eq!(a.transpose_apply(&[1.0, 1.0, 1.0]), vec![4.0, 3.0]);
    }

    /// The aligned backing store must not disturb the row stride.
    ///
    /// A row here is 7 `f64` = 56 bytes, so it neither divides nor is divided by
    /// the 64-byte alignment: consecutive rows straddle cache lines, and 37 rows
    /// span many of them. Both kernels are compared against a naive row-major
    /// reference that sums in the same order, so the comparison is exact.
    #[test]
    fn aligned_storage_preserves_the_row_stride() {
        let (voxels, beamlets) = (37, 7);
        let data: Vec<f64> = (0..voxels * beamlets)
            .map(|i| (i * 37 % 101) as f64 - 50.0)
            .collect();
        let a = DoseInfluence::from_rows(voxels, beamlets, data.clone()).unwrap();
        assert_eq!(a.rows(), data.as_slice());

        let x: Vec<f64> = (0..beamlets).map(|j| j as f64 - 3.0).collect();
        let expected: Vec<f64> = (0..voxels)
            .map(|i| (0..beamlets).map(|j| data[i * beamlets + j] * x[j]).sum())
            .collect();
        assert_eq!(a.apply(&x), expected);

        let r: Vec<f64> = (0..voxels).map(|i| i as f64 * 0.5 - 1.0).collect();
        let expected_transpose: Vec<f64> = (0..beamlets)
            .map(|j| (0..voxels).map(|i| data[i * beamlets + j] * r[i]).sum())
            .collect();
        assert_eq!(a.transpose_apply(&r), expected_transpose);
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

    /// Distance from the exact solution tolerated after a fixed iteration budget.
    ///
    /// Unlike the rounding tolerances elsewhere in this workspace, this one is
    /// **not** derived from `T::EPSILON`: the residual here is algorithmic, not
    /// numerical. With an identity influence matrix each weight converges
    /// geometrically as `(1 - step)^n`, so 300 steps at `step = 0.5` leave a
    /// mathematical residual far below any float width, and the same bound is
    /// therefore correct for every shipped width. Deriving it from epsilon would
    /// assert a precision the iteration count does not deliver, and tightening it
    /// toward epsilon would test the step schedule rather than the optimizer.
    const CONVERGENCE_TOLERANCE: f64 = 1.0e-3;

    /// Asserts the optimizer recovers a separable prescription in one width.
    fn optimizer_recovers_a_separable_prescription<T: ShippedScalar + ScratchElement>() {
        let zero = T::from_f64(0.0);
        let one = T::from_f64(1.0);

        // Identity influence: each beamlet drives exactly one voxel, so the
        // optimum is the prescription itself.
        let mut rows = vec![zero; 4];
        rows[0] = one;
        rows[3] = one;
        let influence =
            DoseInfluence::from_rows(2, 2, rows).expect("2x2 influence matrix is well formed");

        let prescription = [T::from_f64(1.5), T::from_f64(0.5)];
        let weights = optimize_beam_weights(&influence, &prescription, 300, T::from_f64(0.5));

        let tolerance = T::from_f64(CONVERGENCE_TOLERANCE);
        for (got, want) in weights.iter().zip(&prescription) {
            assert_relative_eq!(got, want, epsilon = tolerance);
        }
    }

    #[test]
    fn optimizer_recovers_a_separable_prescription_in_single_precision() {
        optimizer_recovers_a_separable_prescription::<f32>();
    }

    #[test]
    fn optimizer_recovers_a_separable_prescription_in_double_precision() {
        optimizer_recovers_a_separable_prescription::<f64>();
    }
}
