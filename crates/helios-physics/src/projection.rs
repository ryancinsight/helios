//! Ray line-integral reduction for photon projection and dose ray-tracing.
//!
//! A forward projector (MVCT) or a dose ray-trace decomposes into two stages:
//! (1) a **geometry** stage that intersects a ray with the voxel grid and emits
//! the sequence of `(attenuation, segment-length)` pairs along the ray, and
//! (2) a **physics reduction** stage that accumulates those segments into an
//! optical depth and a transmitted fraction.
//!
//! This module owns stage (2) — the reduction — which is pure, geometry-free,
//! and analytically verifiable. Stage (1) (voxel DDA / Siddon traversal) depends
//! on the gaia geometry kernel's `Ray`/`Aabb` and lands once that is available;
//! keeping the reduction separate lets the physics be tested exactly now and lets
//! the same reduction run over CPU- or GPU-generated segments unchanged.
//!
//! The optical depth is the line integral of the linear attenuation coefficient
//! along the ray,
//!
//! ```text
//! τ = ∫ μ(s) ds ≈ Σᵢ μᵢ · Lᵢ ,
//! ```
//!
//! and the narrow-beam transmitted fraction is `exp(−τ)` (Beer–Lambert composed
//! over the path).

use crate::attenuation::LinearAttenuation;
use helios_math::{NumericElement, Scalar};

/// Accumulate optical depth `τ = Σ μᵢ·Lᵢ` over path segments.
///
/// Each segment is `(μ, length_cm)` with `μ` in cm⁻¹ and `length_cm` the
/// geometric path through that segment in cm (expected non-negative). An empty
/// path has zero optical depth.
pub fn optical_depth<T, I>(segments: I) -> T
where
    T: Scalar,
    I: IntoIterator<Item = (LinearAttenuation<T>, T)>,
{
    segments
        .into_iter()
        .fold(<T as NumericElement>::ZERO, |acc, (mu, length_cm)| {
            acc + mu.get() * length_cm
        })
}

/// Narrow-beam transmitted fraction `exp(−τ)` along a path of `(μ, length)`
/// segments.
///
/// Composes the per-segment Beer–Lambert factors: `Π exp(−μᵢ·Lᵢ) = exp(−Σ μᵢ·Lᵢ)`.
/// An empty path transmits fully (`1`).
pub fn beam_transmission<T, I>(segments: I) -> T
where
    T: Scalar,
    I: IntoIterator<Item = (LinearAttenuation<T>, T)>,
{
    (-optical_depth(segments)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn mu(v: f64) -> LinearAttenuation<f64> {
        LinearAttenuation::new(v).expect("valid coefficient")
    }

    #[test]
    fn empty_path_has_zero_depth_and_full_transmission() {
        let empty: [(LinearAttenuation<f64>, f64); 0] = [];
        assert_eq!(optical_depth(empty), 0.0);
        let empty2: [(LinearAttenuation<f64>, f64); 0] = [];
        assert_eq!(beam_transmission(empty2), 1.0);
    }

    #[test]
    fn homogeneous_path_equals_single_slab_mu_l() {
        // Splitting a homogeneous slab into N equal segments must give τ = μ·L
        // exactly (up to summation rounding) — the discretization oracle.
        let mu_val = 0.12;
        let total_len = 7.5;
        let n = 15;
        let seg = total_len / n as f64;
        let segments: Vec<_> = (0..n).map(|_| (mu(mu_val), seg)).collect();
        assert_relative_eq!(optical_depth(segments), mu_val * total_len, epsilon = 1e-12);

        let segments: Vec<_> = (0..n).map(|_| (mu(mu_val), seg)).collect();
        assert_relative_eq!(
            beam_transmission(segments),
            (-(mu_val * total_len)).exp(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn optical_depth_is_additive_over_heterogeneous_segments() {
        // τ over a bone/water/air-like path equals the sum of segment depths.
        let segments = [(mu(0.5), 2.0), (mu(0.18), 3.0), (mu(0.001), 10.0)];
        let expected = 0.5 * 2.0 + 0.18 * 3.0 + 0.001 * 10.0;
        assert_relative_eq!(optical_depth(segments), expected, epsilon = 1e-12);
    }

    #[test]
    fn transmission_composes_multiplicatively() {
        // Transmission of the concatenated path equals the product of the parts.
        let part_a = [(mu(0.3), 4.0)];
        let part_b = [(mu(0.07), 6.0)];
        let whole = [(mu(0.3), 4.0), (mu(0.07), 6.0)];
        assert_relative_eq!(
            beam_transmission(whole),
            beam_transmission(part_a) * beam_transmission(part_b),
            epsilon = 1e-12
        );
    }

    #[test]
    fn reduction_is_generic_over_scalar_f32() {
        let segments = [
            (LinearAttenuation::new(0.2_f32).unwrap(), 3.0_f32),
            (LinearAttenuation::new(0.1_f32).unwrap(), 2.0_f32),
        ];
        // τ = 0.2·3 + 0.1·2 = 0.8.
        assert_relative_eq!(optical_depth(segments), 0.8_f32, epsilon = 1e-6);
    }
}
