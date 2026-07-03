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
//! gradient `Aᵀ(A x − d)` is exact here; a `coeus`-autodiff backend generalizes
//! it to non-quadratic (DVH/biological) objectives in a later increment.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

#[cfg(feature = "autodiff")]
mod autodiff;
mod optimize;

#[cfg(feature = "autodiff")]
pub use autodiff::objective_gradient_autodiff;
pub use optimize::{objective_value, optimize_beam_weights, DoseInfluence};
