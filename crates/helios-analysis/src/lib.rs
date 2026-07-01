//! Helios dosimetric analysis.
//!
//! Plan-quality and plan-comparison metrics over dose [`Volume`](helios_domain::Volume)s:
//! the cumulative [`Dvh`] (dose-volume histogram) and the [`gamma_index_3d`]
//! gamma comparison (Low's combined dose-difference / distance-to-agreement
//! criterion). These are the quantitative gates used to validate dose engines
//! against reference Monte Carlo / measurement (e.g. 3%/2 mm gamma, DVH metrics).
//!
//! All metrics are generic over the [`Scalar`](helios_math::Scalar) seam and
//! operate on the CPU; they are authored independently of the dose engines so the
//! validation machinery is ready as those engines land.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod dvh;
mod gamma;

pub use dvh::Dvh;
pub use gamma::{gamma_index_3d, gamma_pass_rate};
