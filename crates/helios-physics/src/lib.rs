//! Helios radiation-interaction physics.
//!
//! Helios owns modality-specific CT calibration and Compton interaction models.
//! Reusable photon attenuation, optical depth, Beer-Lambert transmission, NIST
//! reference data, and derived optical transport laws are owned by `hyperion`.
//!
//! All quantities are generic over the [`Scalar`](helios_math::Scalar) seam, so a
//! kernel runs natively at `f32` (GPU staging) or `f64` (reference) with no
//! widen-narrow.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

pub mod compton;
mod ct_calibration;

pub use compton::{
    compton_energy_transfer_cross_section, compton_mass_attenuation, compton_mass_energy_transfer,
    compton_mean_energy_transfer_fraction, electrons_per_gram, klein_nishina_cross_section,
    klein_nishina_differential, thomson_cross_section,
};
pub use ct_calibration::{mass_density_from_hu, relative_electron_density_from_hu};
