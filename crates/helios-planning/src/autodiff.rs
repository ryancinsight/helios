//! coeus-autograd gradient backend for the planning objective.
//!
//! Computes `∇ₓ ½‖A·x − d‖²` by reverse-mode automatic differentiation over the
//! coeus tape (the mandated Atlas tensor/autodiff component) instead of the
//! hand-derived `Aᵀ(A·x − d)` in [`optimize_beam_weights`]. For the quadratic objective the two are
//! mathematically identical — the differential test pins them against each other
//! — and the autodiff path is what generalizes to non-quadratic (DVH/biological)
//! objectives, where no closed-form gradient exists.
//!
//! Feature-gated behind `autodiff` so the tensor/tape machinery stays out of the
//! core build. The feature gates a complete implementation, not a stub.

use crate::optimize::DoseInfluence;
use aequitas::systems::si::quantities::{AbsorbedDose, Dimensionless};
use asclepius::VolumeEffect;
use asclepius_coeus::response::radiation::generalized_equivalent_uniform_dose;
use coeus_autograd::{add, matmul, mul, relu, sub, sum, Var};
use coeus_core::MoiraiBackend;
use coeus_tensor::Tensor;
use helios_core::HeliosError;

type AutodiffVar = Var<f64, MoiraiBackend>;

/// A constant (non-differentiated) scalar `Var` of shape `[1]` on `backend`.
#[expect(
    clippy::trivially_copy_pass_by_ref,
    reason = "Tensor::from_slice_on borrows the backend; callers keep ownership"
)]
fn scalar_const(value: f64, backend: &MoiraiBackend) -> AutodiffVar {
    Var::new(Tensor::from_slice_on(vec![1], &[value], backend), false)
}

#[inline]
fn require_beamlet_weights(
    x: &[f64],
    beamlets: usize,
    field: &'static str,
) -> Result<(), HeliosError> {
    if x.len() == beamlets {
        return Ok(());
    }
    Err(HeliosError::InvalidDomainValue {
        field,
        value: x.len() as f64,
        reason: "weight count must equal the beamlet count",
    })
}

#[inline]
fn design_and_weight_var(
    influence: &DoseInfluence<f64>,
    x: &[f64],
    backend: &MoiraiBackend,
) -> (AutodiffVar, AutodiffVar) {
    let (voxels, beamlets) = influence.dims();
    let a = Var::new(
        Tensor::from_slice_on(vec![voxels, beamlets], influence.rows(), backend),
        false,
    );
    let xv = Var::new(Tensor::from_slice_on(vec![beamlets, 1], x, backend), true);
    (a, xv)
}

#[inline]
fn require_gradient(xv: &AutodiffVar, field: &'static str) -> Result<Vec<f64>, HeliosError> {
    let grad = xv.grad().ok_or(HeliosError::InvalidDomainValue {
        field,
        value: f64::NAN,
        reason: "autograd tape produced no gradient for x",
    })?;
    Ok(grad.as_slice().to_vec())
}

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
    require_beamlet_weights(x, beamlets, "objective_gradient_autodiff::x")?;
    if prescription.len() != voxels {
        return Err(HeliosError::InvalidDomainValue {
            field: "objective_gradient_autodiff::prescription",
            value: prescription.len() as f64,
            reason: "prescription length must equal the voxel count",
        });
    }

    let backend = MoiraiBackend::new();
    // A: constants (no gradient tracked); x: the differentiated variable.
    let (a, xv) = design_and_weight_var(influence, x, &backend);
    let d = Var::new(
        Tensor::from_slice_on(vec![voxels, 1], prescription, &backend),
        false,
    );

    let r = sub(&matmul(&a, &xv), &d);
    let loss = sum(&mul(&r, &r)); // ‖r‖² — tape gradient wrt x is 2·Aᵀr.
    let _ = loss.backward();

    let grad = require_gradient(&xv, "objective_gradient_autodiff::grad")?;
    Ok(grad.into_iter().map(|g| 0.5 * g).collect())
}

/// One-sided DVH-style penalty band for the non-quadratic planning objective:
/// underdose below `floor` and overdose above `ceiling` are penalized
/// quadratically; dose inside the band costs nothing.
#[derive(Debug, Clone, Copy)]
pub struct DvhPenalty<'a> {
    /// Per-voxel prescription floor (underdose below it is penalized).
    pub floor: &'a [AbsorbedDose<f64>],
    /// Per-voxel dose ceiling (overdose above it is penalized).
    pub ceiling: &'a [AbsorbedDose<f64>],
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
    if x.len() != beamlets {
        return Err(HeliosError::InvalidDomainValue {
            field: "dvh_objective_gradient_autodiff",
            value: x.len() as f64,
            reason: "weight count must equal the beamlet count",
        });
    }
    for len in [penalty.floor.len(), penalty.ceiling.len()] {
        if len != voxels {
            return Err(HeliosError::InvalidDomainValue {
                field: "dvh_objective_gradient_autodiff",
                value: len as f64,
                reason: "penalty band length must equal the voxel count",
            });
        }
    }

    let backend = MoiraiBackend::new();
    let (a, xv) = design_and_weight_var(influence, x, &backend);
    // The autodiff tensor is a scalar numerical boundary. Preserve the
    // physical dose type through the public contract and unwrap it once here.
    let floor_values: Vec<f64> = penalty.floor.iter().map(|dose| *dose.as_base()).collect();
    let ceiling_values: Vec<f64> = penalty.ceiling.iter().map(|dose| *dose.as_base()).collect();
    let floor = Var::new(
        Tensor::from_slice_on(vec![voxels, 1], &floor_values, &backend),
        false,
    );
    let ceiling = Var::new(
        Tensor::from_slice_on(vec![voxels, 1], &ceiling_values, &backend),
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
    let _ = loss.backward();

    let grad = require_gradient(&xv, "dvh_objective_gradient_autodiff::grad")?;
    Ok((value, grad))
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

/// Which side of the gEUD reference an [`EudPenalty`] penalizes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EudKind {
    /// Penalize gEUD **above** the reference (an OAR dose ceiling; `a ≥ 1`).
    UpperLimit,
    /// Penalize gEUD **below** the reference (a target dose floor; `a < 1`).
    LowerLimit,
}

/// A one-sided quadratic gEUD objective for one structure: penalize the
/// structure's gEUD for violating `reference` on the [`kind`](Self::kind) side,
/// weighted by `weight`.
#[derive(Debug, Clone, Copy)]
pub struct EudPenalty {
    /// Niemierko volume-effect parameter (`a ≠ 0`).
    pub a: Dimensionless<f64>,
    /// gEUD reference dose (ceiling for `UpperLimit`, floor for `LowerLimit`).
    pub reference: AbsorbedDose<f64>,
    /// Which side is penalized.
    pub kind: EudKind,
    /// Penalty weight.
    pub weight: f64,
}

/// Objective value and gradient of a one-sided gEUD penalty
/// `weight · relu(±(gEUD(A·x) − reference))²` with respect to the beam weights
/// `x`, computed on the coeus tape.
///
/// The dose field is passed directly to the differentiable Asclepius gEUD law,
/// so Helios owns only the planning objective while Asclepius owns its response
/// equation and stabilized tape construction. The tests pin the returned
/// gradient against a central finite-difference of the objective.
///
/// # Errors
/// [`HeliosError::InvalidDomainValue`] if `x`'s length differs from the beamlet
/// count, or `penalty.a == 0`.
pub fn eud_objective_gradient_autodiff(
    influence: &DoseInfluence<f64>,
    x: &[f64],
    penalty: &EudPenalty,
) -> Result<(f64, Vec<f64>), HeliosError> {
    let (_, beamlets) = influence.dims();
    require_beamlet_weights(x, beamlets, "eud_objective_gradient_autodiff::x")?;
    let volume_effect =
        VolumeEffect::new(*penalty.a.as_base()).map_err(|_| HeliosError::InvalidDomainValue {
            field: "eud_objective_gradient_autodiff::a",
            value: *penalty.a.as_base(),
            reason: "gEUD volume parameter must be finite and non-zero",
        })?;

    let backend = MoiraiBackend::new();
    let (a, xv) = design_and_weight_var(influence, x, &backend);

    let dose = matmul(&a, &xv);
    let geud = generalized_equivalent_uniform_dose(&dose, volume_effect).map_err(|error| {
        let value = match error {
            asclepius_coeus::AutodiffResponseError::InvalidDose { value, .. } => value,
            _ => f64::NAN,
        };
        HeliosError::InvalidDomainValue {
            field: "eud_objective_gradient_autodiff::dose",
            value,
            reason: "dose observation violates the Asclepius gEUD domain",
        }
    })?;

    // One-sided hinge: violation = ±(gEUD − reference), penalized when > 0.
    let reference = scalar_const(*penalty.reference.as_base(), &backend);
    let violation = match penalty.kind {
        EudKind::UpperLimit => sub(&geud, &reference),
        EudKind::LowerLimit => sub(&reference, &geud),
    };
    let hinge = relu(&violation);
    let weight = scalar_const(penalty.weight, &backend);
    let loss = mul(&mul(&hinge, &hinge), &weight);

    let value = loss.tensor.as_slice()[0];
    let _ = loss.backward();
    let grad = require_gradient(&xv, "eud_objective_gradient_autodiff::grad")?;
    Ok((value, grad))
}

#[cfg(test)]
mod tests;
