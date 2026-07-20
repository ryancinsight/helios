//! MVCT quantum (photon-counting) noise model for projection sinograms.
//!
//! A physical MVCT detector measures a transmitted photon count
//! `N = N₀·exp(−τ)` for an incident fluence of `N₀` photons per ray and a line
//! integral `τ`. Photon counting is Poisson, so `N` carries variance `N`; the
//! reconstructed line integral `τ' = −ln(N'/N₀)` inherits noise that **grows with
//! attenuation** — the hallmark of CT/MVCT quantum noise. This module injects that
//! noise deterministically (seeded) so noisy reconstructions can be produced and
//! the [`image_quality`](../../helios_analysis/image_quality/index.html) noise/CNR
//! metrics exercised end-to-end.
//!
//! Photon counts in the diagnostic/MV regime are large (`N ≫ 20`), so the Poisson
//! draw uses the Gaussian approximation `N' = N + √N·z`, `z ~ 𝒩(0,1)` (exact in
//! the large-count limit; error `O(1/√N)`). Sampling is done in `f64` — the noise
//! model's natural precision — and cast back to the sinogram's scalar type.

use crate::radon::Sinogram;
use helios_math::GeometryScalar;
use tyche_core::{Seed, StandardNormal};

/// Add MVCT quantum noise to a sinogram of line integrals `τ`, returning the
/// noise-perturbed sinogram.
///
/// For each reading: `N = photons_per_ray·exp(−τ)`, draw `N' = N + √N·z`
/// (`z ~ 𝒩(0,1)`), clamp to `≥ 1` photon, and return `τ' = −ln(N'/photons_per_ray)`.
/// By error propagation the per-reading noise variance is `Var(τ') ≈ exp(τ)/N₀`,
/// so noise rises with attenuation and falls as `1/N₀` — the analytical oracle in
/// the tests. Deterministic in `seed`; row-major `[angle][offset]` order maps
/// directly to Tyche's logical sample index.
///
/// # Panics
/// Does not panic; `photons_per_ray` should be positive (a non-positive value
/// yields a degenerate all-`NaN`/`Inf` sinogram, the caller's contract).
#[must_use]
pub fn add_quantum_noise<T: GeometryScalar>(
    sinogram: &Sinogram<T>,
    photons_per_ray: f64,
    seed: u64,
) -> Sinogram<T> {
    let seed = Seed::new(seed);
    let mut sample_index = 0_u64;
    sinogram.map_readings(|tau_t| {
        let tau = tau_t.to_f64();
        let expected = photons_per_ray * (-tau).exp();
        let gaussian = StandardNormal::at(seed, sample_index, 0);
        sample_index = sample_index.wrapping_add(1);
        let noisy = expected + expected.sqrt() * gaussian;
        let counts = noisy.max(1.0);
        <T as GeometryScalar>::from_f64(-(counts / photons_per_ray).ln())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Constant-τ sinogram of `n` readings (1 angle × n offsets), for statistics.
    fn constant_sinogram(tau: f64, n: usize) -> Sinogram<f64> {
        let offsets: Vec<f64> = (0..n).map(|i| i as f64).collect();
        Sinogram::from_readings(vec![0.0], offsets, vec![tau; n]).expect("valid sinogram")
    }

    fn mean_and_var(sino: &Sinogram<f64>) -> (f64, f64) {
        let (_, n_off) = sino.dims();
        let vals: Vec<f64> = (0..n_off).map(|j| sino.get(0, j)).collect();
        let n = vals.len() as f64;
        let mean = vals.iter().sum::<f64>() / n;
        let var = vals.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
        (mean, var)
    }

    #[test]
    fn same_seed_is_deterministic_different_seed_differs() {
        let clean = constant_sinogram(0.5, 256);
        let a = add_quantum_noise(&clean, 1.0e4, 42);
        let b = add_quantum_noise(&clean, 1.0e4, 42);
        let c = add_quantum_noise(&clean, 1.0e4, 43);
        assert_eq!(a, b, "same seed must reproduce exactly");
        assert_ne!(a, c, "different seed must differ");
    }

    #[test]
    fn tyche_seed_mapping_is_pinned() {
        let clean = constant_sinogram(0.5, 1);
        let noisy = add_quantum_noise(&clean, 1.0e4, 42);
        assert_eq!(noisy.get(0, 0).to_bits(), 0x3FDF_467D_FB82_DEC9);
    }

    #[test]
    fn high_flux_limit_recovers_clean_line_integrals() {
        // As N₀ → ∞ the relative noise 1/√N → 0, so τ' → τ.
        let clean = constant_sinogram(0.7, 128);
        let noisy = add_quantum_noise(&clean, 1.0e14, 7);
        for j in 0..128 {
            assert!(
                (noisy.get(0, j) - 0.7).abs() < 1e-5,
                "reading {j} not ~clean"
            );
        }
    }

    #[test]
    fn noise_variance_matches_photon_statistics() {
        // Var(τ') ≈ exp(τ)/N₀. Ensemble of n independent readings; sample variance
        // has relative error √(2/n) ≈ 1% at n = 20000 → 10% tolerance is safe.
        let tau0 = 0.5;
        let n0 = 1.0e4;
        let noisy = add_quantum_noise(&constant_sinogram(tau0, 20_000), n0, 2024);
        let (mean, var) = mean_and_var(&noisy);
        let expected_var = tau0.exp() / n0;
        assert!(
            (var / expected_var - 1.0).abs() < 0.10,
            "sample var {var:.3e} vs analytical {expected_var:.3e}"
        );
        // Mean is ~unbiased (bias ~ 1/2N ≪ std/√n).
        assert!((mean - tau0).abs() < 3e-3, "mean {mean} biased from {tau0}");
    }

    #[test]
    fn noise_grows_with_attenuation() {
        // Var(τ') = exp(τ)/N₀ increases with τ: a thicker path is noisier.
        let n0 = 1.0e4;
        let (_, var_low) = mean_and_var(&add_quantum_noise(&constant_sinogram(0.2, 20_000), n0, 1));
        let (_, var_high) =
            mean_and_var(&add_quantum_noise(&constant_sinogram(1.5, 20_000), n0, 1));
        assert!(
            var_high > var_low * 2.0,
            "high-τ noise {var_high:.3e} should dominate low-τ {var_low:.3e}"
        );
    }

    #[test]
    fn noise_model_is_generic_over_scalar_f32() {
        let offsets: Vec<f32> = (0..64).map(|i| i as f32).collect();
        let clean = Sinogram::from_readings(vec![0.0_f32], offsets, vec![0.5_f32; 64]).unwrap();
        let a = add_quantum_noise(&clean, 1.0e4, 99);
        let b = add_quantum_noise(&clean, 1.0e4, 99);
        assert_eq!(a, b);
        // High flux → near-clean.
        let hi = add_quantum_noise(&clean, 1.0e12, 99);
        assert!((hi.get(0, 0) - 0.5_f32).abs() < 1e-3);
    }
}
