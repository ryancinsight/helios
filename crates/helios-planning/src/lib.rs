//! Helios inverse treatment planning.
//!
//! Fluence-map / beam-weight optimization: given a linear dose model
//! `dose = A · x` (the [`DoseInfluence`] matrix mapping beamlet weights `x` to
//! voxel dose) and a prescribed dose `d`, find non-negative weights minimizing
//! the quadratic objective `½‖A x − d‖²` by projected gradient descent
//! ([`optimize_beam_weights`]).
//!
//! The objective is convex, so the projected gradient converges to the
//! constrained least-squares optimum — the analytical oracle in the tests. The
//! gradient `Aᵀ(A x − d)` is exact here; the `coeus`-autodiff backend (feature
//! `autodiff`) generalizes it to non-quadratic objectives with no closed-form
//! gradient — DVH-band penalties and generalized-EUD (Niemierko) biological
//! objectives.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

#[cfg(feature = "autodiff")]
mod autodiff;
mod optimize;

#[cfg(feature = "autodiff")]
pub use autodiff::{
    dvh_objective_gradient_autodiff, eud_objective_gradient_autodiff, objective_gradient_autodiff,
    optimize_beam_weights_dvh, DvhPenalty, EudKind, EudPenalty,
};
pub use optimize::{objective_value, optimize_beam_weights, DoseInfluence};
