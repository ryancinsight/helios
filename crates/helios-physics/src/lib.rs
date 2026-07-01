//! Helios radiation-interaction physics.
//!
//! Photon transport primitives shared by imaging (MVCT forward projection,
//! reconstruction) and therapy (dose ray-tracing): the [`attenuation`] module
//! provides validated linear/mass attenuation coefficients, the Beer–Lambert
//! transmission law, and CT-number-to-density calibration.
//!
//! All quantities are generic over the [`Scalar`](helios_math::Scalar) seam, so a
//! kernel runs natively at `f32` (GPU staging) or `f64` (reference) with no
//! widen-narrow. Concrete cross-section tables (e.g. NIST XCOM μ/ρ) are data that
//! later increments load; this crate owns the physics relations they feed.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

pub mod attenuation;
pub mod compton;
pub mod projection;

pub use attenuation::{
    mass_density_from_hu, relative_electron_density_from_hu, LinearAttenuation, MassAttenuation,
};
pub use compton::{
    compton_energy_transfer_cross_section, compton_mass_attenuation, compton_mass_energy_transfer,
    compton_mean_energy_transfer_fraction, electrons_per_gram, klein_nishina_cross_section,
    klein_nishina_differential, thomson_cross_section,
};
pub use projection::{beam_transmission, optical_depth};
