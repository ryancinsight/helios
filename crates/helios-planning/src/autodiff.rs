//! coeus-autograd gradient backend for the planning objective.
//!
//! Computes `∇ₓ ½‖A·x − d‖²` by reverse-mode automatic differentiation over the
//! coeus tape (the mandated Atlas tensor/autodiff component) instead of the
//! hand-derived `Aᵀ(A·x − d)` in [`optimize_beam_weights`]
//! (crate::optimize_beam_weights). For the quadratic objective the two are
//! mathematically identical — the differential test pins them against each other
//! — and the autodiff path is what generalizes to non-quadratic (DVH/biological)
//! objectives, where no closed-form gradient exists.
//!
//! Feature-gated behind `autodiff` so the tensor/tape machinery stays out of the
//! core build. The feature gates a complete implementation, not a stub.

use crate::optimize::DoseInfluence;
use coeus_autograd::{matmul, mul, sub, sum, Var};
use coeus_core::MoiraiBackend;
use coeus_tensor::Tensor;
use helios_core::HeliosError;

/// Gradient of the quadratic objective `½‖A·x − d‖²` with respect to `x`,
/// computed by coeus reverse-mode autodiff.
///
/// Tape: `r = A·x − d` → `loss = Σ r⊙r = ‖r‖²` → backward. The tape gradient of
/// `‖r‖²` is `2·Aᵀr`, so the returned gradient is halved to match the
/// `½‖·‖²` convention (and the exact hand gradient in
/// [`optimize_beam_weights`](crate::optimize_beam_weights)).
///
/// Values cross the coeus boundary as `f64` (the tensor backend's reference
/// precision; the same concrete-at-the-boundary convention as the PyO3/DICOM
/// boundaries).
///
/// # Errors
/// [`HeliosError::InvalidDomainValue`] if `x`/`prescription` lengths do not match
/// `influence.dims()`.
pub fn objective_gradient_autodiff(
    influence: &DoseInfluence<f64>,
    x: &[f64],
    prescription: &[f64],
) -> Result<Vec<f64>, HeliosError> {
    let (voxels, beamlets) = influence.dims();
    if x.len() != beamlets {
        return Err(HeliosError::InvalidDomainValue {
            field: "objective_gradient_autodiff::x",
            value: x.len() as f64,
            reason: "weight count must equal the beamlet count",
        });
    }
    if prescription.len() != voxels {
        return Err(HeliosError::InvalidDomainValue {
            field: "objective_gradient_autodiff::prescription",
            value: prescription.len() as f64,
            reason: "prescription length must equal the voxel count",
        });
    }

    let backend = MoiraiBackend::new();
    // A: constants (no gradient tracked); x: the differentiated variable.
    let a = Var::<f64, MoiraiBackend>::new(
        Tensor::from_slice_on(vec![voxels, beamlets], influence.rows(), &backend),
        false,
    );
    let xv = Var::new(Tensor::from_slice_on(vec![beamlets, 1], x, &backend), true);
    let d = Var::new(
        Tensor::from_slice_on(vec![voxels, 1], prescription, &backend),
        false,
    );

    let r = sub(&matmul(&a, &xv), &d);
    let loss = sum(&mul(&r, &r)); // ‖r‖² — tape gradient wrt x is 2·Aᵀr.
    loss.backward();

    let grad = xv.grad().ok_or(HeliosError::InvalidDomainValue {
        field: "objective_gradient_autodiff::grad",
        value: f64::NAN,
        reason: "autograd tape produced no gradient for x",
    })?;
    Ok(grad.as_slice().iter().map(|&g| 0.5 * g).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// 3×2 influence with distinct entries; the differential oracle is the exact
    /// hand gradient Aᵀ(A·x − d) from the projected-gradient solver.
    fn influence() -> DoseInfluence<f64> {
        DoseInfluence::from_rows(3, 2, vec![1.0, 2.0, 0.5, -1.0, 3.0, 4.0]).unwrap()
    }

    #[test]
    fn autodiff_gradient_matches_the_exact_hand_gradient() {
        let inf = influence();
        let x = [0.7, -0.3];
        let d = [1.0, 0.5, 2.0];

        // Exact: Aᵀ(A·x − d).
        let ax = inf.apply(&x);
        let residual: Vec<f64> = ax.iter().zip(&d).map(|(&a, &b)| a - b).collect();
        let exact = inf.transpose_apply(&residual);

        let auto = objective_gradient_autodiff(&inf, &x, &d).expect("gradient");
        assert_eq!(auto.len(), exact.len());
        for (i, (&g, &e)) in auto.iter().zip(&exact).enumerate() {
            // Same values through a different summation route: bound at 1e-12
            // relative (f64; tiny fixed dims, reduction depth ≤ 3).
            assert_relative_eq!(g, e, max_relative = 1e-12, epsilon = 1e-14);
            let _ = i;
        }
    }

    #[test]
    fn gradient_is_zero_at_the_least_squares_optimum() {
        // Identity influence: optimum x* = d exactly → ∇ = 0.
        let inf = DoseInfluence::from_rows(2, 2, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
        let d = [2.0, -1.5];
        let grad = objective_gradient_autodiff(&inf, &d, &d).expect("gradient");
        for &g in &grad {
            assert_relative_eq!(g, 0.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn shape_mismatches_are_typed_errors() {
        let inf = influence();
        assert!(objective_gradient_autodiff(&inf, &[1.0], &[1.0, 1.0, 1.0]).is_err());
        assert!(objective_gradient_autodiff(&inf, &[1.0, 2.0], &[1.0]).is_err());
    }
}
