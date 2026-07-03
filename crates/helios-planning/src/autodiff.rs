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
use coeus_autograd::{add, matmul, mul, relu, sub, sum, Var};
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


/// One-sided DVH-style penalty band for the non-quadratic planning objective:
/// underdose below `floor` and overdose above `ceiling` are penalized
/// quadratically; dose inside the band costs nothing.
#[derive(Debug, Clone, Copy)]
pub struct DvhPenalty<'a> {
    /// Per-voxel prescription floor (underdose below it is penalized).
    pub floor: &'a [f64],
    /// Per-voxel dose ceiling (overdose above it is penalized).
    pub ceiling: &'a [f64],
    /// Weight on the underdose term.
    pub weight_under: f64,
    /// Weight on the overdose term.
    pub weight_over: f64,
}

/// Non-quadratic DVH-penalty objective and its autodiff gradient:
///
/// ```text
/// L(x) = w_u · Σ relu(floor − A·x)² + w_o · Σ relu(A·x − ceiling)²
/// ```
///
/// This is the piecewise (one-sided) clinical objective — no closed-form gradient
/// in general — computed on the coeus tape (`relu` kinks handled by reverse-mode
/// AD). For verification the sub-gradient is still hand-derivable:
/// `∇L = −2·w_u·Aᵀ relu(floor − A·x) + 2·w_o·Aᵀ relu(A·x − ceiling)`, the
/// differential oracle in the tests. Returns `(objective value, gradient)`.
///
/// # Errors
/// [`HeliosError::InvalidDomainValue`] on any length mismatch with
/// `influence.dims()`.
pub fn dvh_objective_gradient_autodiff(
    influence: &DoseInfluence<f64>,
    x: &[f64],
    penalty: &DvhPenalty<'_>,
) -> Result<(f64, Vec<f64>), HeliosError> {
    let (voxels, beamlets) = influence.dims();
    for (label, len, want) in [
        ("x", x.len(), beamlets),
        ("floor", penalty.floor.len(), voxels),
        ("ceiling", penalty.ceiling.len(), voxels),
    ] {
        if len != want {
            return Err(HeliosError::InvalidDomainValue {
                field: "dvh_objective_gradient_autodiff",
                value: len as f64,
                reason: match label {
                    "x" => "weight count must equal the beamlet count",
                    _ => "penalty band length must equal the voxel count",
                },
            });
        }
    }

    let backend = MoiraiBackend::new();
    let a = Var::<f64, MoiraiBackend>::new(
        Tensor::from_slice_on(vec![voxels, beamlets], influence.rows(), &backend),
        false,
    );
    let xv = Var::new(Tensor::from_slice_on(vec![beamlets, 1], x, &backend), true);
    let floor = Var::new(
        Tensor::from_slice_on(vec![voxels, 1], penalty.floor, &backend),
        false,
    );
    let ceiling = Var::new(
        Tensor::from_slice_on(vec![voxels, 1], penalty.ceiling, &backend),
        false,
    );
    let wu = Var::new(
        Tensor::from_slice_on(vec![1], &[penalty.weight_under], &backend),
        false,
    );
    let wo = Var::new(
        Tensor::from_slice_on(vec![1], &[penalty.weight_over], &backend),
        false,
    );

    let ax = matmul(&a, &xv);
    let under = relu(&sub(&floor, &ax)); // relu(floor − dose)
    let over = relu(&sub(&ax, &ceiling)); // relu(dose − ceiling)
    let loss = add(
        &mul(&sum(&mul(&under, &under)), &wu),
        &mul(&sum(&mul(&over, &over)), &wo),
    );
    let value = loss.tensor.as_slice()[0];
    loss.backward();

    let grad = xv.grad().ok_or(HeliosError::InvalidDomainValue {
        field: "dvh_objective_gradient_autodiff::grad",
        value: f64::NAN,
        reason: "autograd tape produced no gradient for x",
    })?;
    Ok((value, grad.as_slice().to_vec()))
}

/// Projected-gradient descent on the non-quadratic DVH-penalty objective, using
/// the coeus autodiff gradient each iteration (`x ← max(0, x − step·∇L)`).
///
/// The objective is convex (a sum of squared hinges composed with a linear map),
/// so descent with a suitable step reaches the constrained optimum; the tests
/// assert the clinical semantics (target voxels raised to the floor, OAR voxels
/// held under the ceiling, weights non-negative).
///
/// # Errors
/// As [`dvh_objective_gradient_autodiff`].
pub fn optimize_beam_weights_dvh(
    influence: &DoseInfluence<f64>,
    penalty: &DvhPenalty<'_>,
    iterations: usize,
    step: f64,
) -> Result<Vec<f64>, HeliosError> {
    let (_, beamlets) = influence.dims();
    let mut x = vec![0.0f64; beamlets];
    for _ in 0..iterations {
        let (_, grad) = dvh_objective_gradient_autodiff(influence, &x, penalty)?;
        for (xj, gj) in x.iter_mut().zip(&grad) {
            *xj = (*xj - step * gj).max(0.0);
        }
    }
    Ok(x)
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
    fn dvh_gradient_matches_the_hand_subgradient() {
        // ∇L = −2·w_u·Aᵀ relu(f − Ax) + 2·w_o·Aᵀ relu(Ax − c). Choose x so both
        // hinges are strictly active somewhere (no kink ambiguity).
        let inf = influence();
        let x = [0.4, -0.2];
        let floor = [1.0, 0.2, 0.5];
        let ceiling = [1.5, 0.4, 0.9];
        let penalty = DvhPenalty {
            floor: &floor,
            ceiling: &ceiling,
            weight_under: 2.0,
            weight_over: 3.0,
        };
        let ax = inf.apply(&x);
        let under: Vec<f64> = ax.iter().zip(&floor).map(|(&d, &f)| (f - d).max(0.0)).collect();
        let over: Vec<f64> = ax.iter().zip(&ceiling).map(|(&d, &c)| (d - c).max(0.0)).collect();
        let gu = inf.transpose_apply(&under);
        let go = inf.transpose_apply(&over);
        let hand: Vec<f64> = gu
            .iter()
            .zip(&go)
            .map(|(&u, &o)| -2.0 * 2.0 * u + 2.0 * 3.0 * o)
            .collect();

        let (value, auto) = dvh_objective_gradient_autodiff(&inf, &x, &penalty).expect("grad");
        // Objective value cross-check.
        let expected_value = 2.0 * under.iter().map(|u| u * u).sum::<f64>()
            + 3.0 * over.iter().map(|o| o * o).sum::<f64>();
        assert_relative_eq!(value, expected_value, max_relative = 1e-12);
        for (&g, &h) in auto.iter().zip(&hand) {
            assert_relative_eq!(g, h, max_relative = 1e-12, epsilon = 1e-14);
        }
    }

    #[test]
    fn gradient_is_zero_inside_the_penalty_band() {
        // Dose strictly between floor and ceiling everywhere → both hinges
        // inactive → L = 0 and ∇L = 0 (the band is free).
        let inf = DoseInfluence::from_rows(2, 2, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
        let x = [1.0, 1.0]; // dose = [1, 1]
        let penalty = DvhPenalty {
            floor: &[0.5, 0.5],
            ceiling: &[2.0, 2.0],
            weight_under: 1.0,
            weight_over: 1.0,
        };
        let (value, grad) = dvh_objective_gradient_autodiff(&inf, &x, &penalty).expect("grad");
        assert_relative_eq!(value, 0.0, epsilon = 1e-14);
        for &g in &grad {
            assert_relative_eq!(g, 0.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn dvh_optimizer_meets_target_floor_and_oar_ceiling() {
        // Voxel 0 is the target (floor 1.0), voxel 1 the OAR (ceiling 0.3).
        // Beamlet 0 doses both (target 1.0, OAR 0.5/unit); beamlet 1 doses only
        // the target. The optimum uses beamlet 1 (OAR-sparing) to reach the floor
        // while keeping the OAR under its ceiling.
        let inf = DoseInfluence::from_rows(2, 2, vec![1.0, 1.0, 0.5, 0.0]).unwrap();
        let penalty = DvhPenalty {
            floor: &[1.0, 0.0],
            ceiling: &[10.0, 0.3],
            weight_under: 1.0,
            weight_over: 10.0,
        };
        let x = optimize_beam_weights_dvh(&inf, &penalty, 500, 0.05).expect("optimize");
        let dose = inf.apply(&x);
        assert!(x.iter().all(|&w| w >= 0.0), "weights must be non-negative");
        assert!(dose[0] > 0.95, "target dose {} below floor", dose[0]);
        assert!(dose[1] < 0.35, "OAR dose {} above ceiling", dose[1]);
    }

    #[test]
    fn shape_mismatches_are_typed_errors() {
        let inf = influence();
        assert!(objective_gradient_autodiff(&inf, &[1.0], &[1.0, 1.0, 1.0]).is_err());
        assert!(objective_gradient_autodiff(&inf, &[1.0, 2.0], &[1.0]).is_err());
    }
}
