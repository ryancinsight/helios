//! Shared transport-quantity conversions.
//!
//! One definition of the projection→transmission step (`τ ↦ exp(−τ)`) so the
//! helical acquisition, portal-dosimetry, and quantum-noise paths cannot drift
//! apart. This is the lowest crate those three share that also owns Hyperion.

use aequitas::systems::si::quantities::Dimensionless;
use helios_math::Scalar;
use hyperion::{quantity::OpticalDepth, TransportError};

/// Beer–Lambert transmission `exp(−τ)` for an optical depth `τ`.
///
/// Thin wrapper over Hyperion's validated [`OpticalDepth`]: a negative or
/// non-finite `tau` is rejected, and the transmitted fraction is returned as a
/// [`Dimensionless`] quantity.
///
/// # Errors
///
/// Returns [`TransportError::InvalidValue`] if `tau` is negative or non-finite.
pub fn transmission_of<T: Scalar>(
    tau: Dimensionless<T>,
) -> Result<Dimensionless<T>, TransportError<T>> {
    Ok(OpticalDepth::new(tau)?.transmission().into_quantity())
}
