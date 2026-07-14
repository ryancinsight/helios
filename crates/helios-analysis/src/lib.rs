//! Helios dosimetric analysis.
//!
//! Plan-quality and plan-comparison metrics over dose [`Volume`](helios_domain::Volume)s:
//! the cumulative [`Dvh`] (dose-volume histogram) and the [`gamma_index_3d`]
//! gamma comparison (Low's combined dose-difference / distance-to-agreement
//! criterion). These are the quantitative gates used to validate dose engines
//! against reference Monte Carlo / measurement (e.g. 3%/2 mm gamma, DVH metrics).
//! The `image_quality` metrics (reconstruction accuracy, noise, contrast/CNR)
//! are the analogous instruments for the MVCT imaging gate, and the
//! `radiobiology` metrics ([`generalized_eud`], [`tcp_logistic`],
//! [`ntcp_lkb`]) collapse a dose distribution to outcome-probability predictions.
//!
//! All metrics are generic over the [`Scalar`](helios_math::Scalar) seam and
//! operate on the CPU; they are authored independently of the dose engines so the
//! validation machinery is ready as those engines land.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod dvh;
mod gamma;
mod image_quality;
mod radiobiology;
mod roi;

pub use dvh::Dvh;
pub use gamma::{gamma_index_3d, gamma_index_3d_local, gamma_pass_rate};
pub use image_quality::{
    contrast_to_noise_ratio, michelson_contrast, roi_statistics, volume_relative_l2_error,
    volume_rmse, RoiStats,
};
pub use radiobiology::{generalized_eud, ntcp_lkb, tcp_logistic};
pub use roi::{box_mask, spherical_mask};
