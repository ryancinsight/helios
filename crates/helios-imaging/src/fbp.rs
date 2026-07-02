//! Filtered back-projection (FBP) reconstruction of a parallel-beam sinogram.
//!
//! Inverts the [`Sinogram`] Radon transform: each projection is convolved with
//! the Ram-Lak ramp filter, then back-projected across all angles onto the
//! reconstruction grid. Detector offsets and the ramp are handled in **cm** to
//! match the projector's cm line integral, so the reconstruction recovers the
//! linear attenuation `μ` (cm⁻¹) directly.
//!
//! Assumes uniformly-spaced projection angles over `[0, π)` and uniformly-spaced
//! detector offsets — the standard parallel-beam FBP sampling.

use crate::radon::Sinogram;
use helios_core::constants::MM_PER_CM;
use helios_domain::{Volume, VoxelGrid};
use helios_math::{GeometryScalar, NumericElement};

/// Ram-Lak ramp-filter kernel `h[n]` for `n ∈ [−(len−1), len−1]`, sample spacing
/// `ds_cm`. `h[0]=1/(4Δs²)`, `h[odd]=−1/(π²n²Δs²)`, `h[even≠0]=0`. Returned as a
/// `2·len−1` vector indexed by `n + (len−1)`.
fn ram_lak_kernel<T: GeometryScalar>(len: usize, ds_cm: T) -> Vec<T> {
    let zero = <T as NumericElement>::ZERO;
    let inv_ds_sq = (ds_cm * ds_cm).recip();
    let quarter = <T as GeometryScalar>::from_f64(0.25);
    let pi_sq = T::PI * T::PI;
    let mut kernel = vec![zero; 2 * len - 1];
    let base = len as isize - 1;
    for n in -base..=base {
        let value = if n == 0 {
            quarter * inv_ds_sq
        } else if n % 2 != 0 {
            let n_t = <T as GeometryScalar>::from_f64((n * n) as f64);
            -(pi_sq * n_t).recip() * inv_ds_sq
        } else {
            zero
        };
        kernel[(n + base) as usize] = value;
    }
    kernel
}

/// Reconstruct the axial `μ` slice from `sinogram` onto `recon` by filtered
/// back-projection.
///
/// The reconstruction grid's axial centre is taken as the rotation centre (the
/// same centre the forward projection used). Values recover linear attenuation
/// `μ` (cm⁻¹).
#[must_use]
pub fn filtered_back_projection<T: GeometryScalar>(
    sinogram: &Sinogram<T>,
    recon: &VoxelGrid<T>,
) -> Volume<T> {
    let zero = <T as NumericElement>::ZERO;
    let (n_ang, n_off) = sinogram.dims();
    let angles = sinogram.angles();
    let offsets = sinogram.offsets();
    let mm_to_cm = <T as GeometryScalar>::from_f64(MM_PER_CM).recip();

    // Detector spacing in cm (drives the ramp-filter sample spacing).
    let ds_cm = (offsets[1] - offsets[0]) * mm_to_cm;
    let kernel = ram_lak_kernel::<T>(n_off, ds_cm);
    let base = n_off as isize - 1;

    // Ram-Lak-filtered projections: filtered[a][i] = Δs · Σ_k p[a][k]·h[i−k].
    let mut filtered = vec![zero; n_ang * n_off];
    for a in 0..n_ang {
        for i in 0..n_off {
            let mut acc = zero;
            for k in 0..n_off {
                let m = i as isize - k as isize; // in [-(n_off-1), n_off-1]
                acc += sinogram.get(a, k) * kernel[(m + base) as usize];
            }
            filtered[a * n_off + i] = acc * ds_cm;
        }
    }

    // Back-project the ramp-filtered rows, weighted by the angular step Δθ.
    let d_theta = if n_ang > 1 {
        angles[1] - angles[0]
    } else {
        T::PI
    };
    crate::backproject::back_project_rows(angles, offsets, &filtered, recon, d_theta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::radon::parallel_beam_radon;
    use approx::assert_relative_eq;
    use helios_math::Point3;

    fn disk_phantom(mu0: f64, radius_mm: f64) -> Volume<f64> {
        let n = 161;
        let spacing = 0.5;
        let grid = VoxelGrid::axis_aligned(
            [n, n, 1],
            [spacing, spacing, spacing],
            Point3::new(0.0, 0.0, 0.0),
        )
        .unwrap();
        let centre = (n - 1) as f64 * spacing / 2.0;
        Volume::from_shape_fn(grid, move |idx| {
            let dx = idx[0] as f64 * spacing - centre;
            let dy = idx[1] as f64 * spacing - centre;
            if (dx * dx + dy * dy).sqrt() <= radius_mm {
                mu0
            } else {
                0.0
            }
        })
    }

    // Recon grid: 2 mm slice, centre aligned with the phantom centre (40 mm).
    fn recon_grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([41, 41, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0)).unwrap()
    }

    fn uniform_angles(n: usize) -> Vec<f64> {
        (0..n)
            .map(|a| a as f64 * std::f64::consts::PI / n as f64)
            .collect()
    }

    fn uniform_offsets(half_mm: f64, n: usize) -> Vec<f64> {
        let ds = 2.0 * half_mm / (n - 1) as f64;
        (0..n).map(|j| -half_mm + j as f64 * ds).collect()
    }

    #[test]
    fn fbp_round_trip_recovers_disk_attenuation() {
        // disk μ=0.04 cm⁻¹, R=25 mm → radon → FBP recovers μ at the centre and
        // ~0 in the background.
        let mu0 = 0.04;
        let phantom = disk_phantom(mu0, 25.0);
        let angles = uniform_angles(180);
        let offsets = uniform_offsets(45.0, 181);
        let sino = parallel_beam_radon(&phantom, &angles, &offsets, 400.0, 0.25);
        let recon = filtered_back_projection(&sino, &recon_grid());

        // Centre pixel (world 40,40) ≈ μ₀ within FBP tolerance.
        let centre = recon.get(20, 20, 0).unwrap();
        assert_relative_eq!(centre, mu0, max_relative = 0.15);
        // A pixel well inside the disk is also ≈ μ₀ (uniform interior).
        let inner = recon.get(23, 20, 0).unwrap(); // world (46,40), ~6 mm off-centre
        assert_relative_eq!(inner, mu0, max_relative = 0.2);
        // Background well outside the disk is near zero.
        let background = recon.get(2, 2, 0).unwrap(); // world (4,4)
        assert!(
            background.abs() < 0.1 * mu0,
            "background {background} not ~0"
        );
    }

    #[test]
    fn mvct_reconstruction_quality_metrics() {
        // End-to-end MVCT image-quality gate: reconstruct the disk and quantify
        // accuracy (interior mean vs μ₀), background suppression, contrast, and CNR
        // with the helios-analysis metrics.
        use helios_analysis::{contrast_to_noise_ratio, michelson_contrast, roi_statistics};
        let mu0 = 0.04;
        let phantom = disk_phantom(mu0, 25.0);
        let angles = uniform_angles(180);
        let offsets = uniform_offsets(45.0, 181);
        let sino = parallel_beam_radon(&phantom, &angles, &offsets, 400.0, 0.25);
        let recon = filtered_back_projection(&sino, &recon_grid());

        // Interior ROI (well inside the 25 mm disk around centre voxel 20,20):
        // reconstruction accuracy — mean recovers μ₀ within FBP tolerance.
        let interior = roi_statistics(&recon, [18, 18, 0], [23, 23, 1]);
        assert_relative_eq!(interior.mean, mu0, max_relative = 0.15);

        // Background ROI (corner, far outside the disk): near zero.
        let background = roi_statistics(&recon, [0, 0, 0], [5, 5, 1]);
        assert!(
            background.mean.abs() < 0.1 * mu0,
            "background mean {} not ~0",
            background.mean
        );

        // Contrast disk-vs-air is near 1; CNR shows the disk is clearly detectable
        // above the interior reconstruction ripple.
        let contrast = michelson_contrast(interior.mean, background.mean.abs());
        assert!(
            contrast > 0.85,
            "disk/background contrast {contrast} too low"
        );
        let cnr = contrast_to_noise_ratio(interior.mean, background.mean, interior.std);
        assert!(cnr > 1.0, "cnr {cnr} indicates the disk is not detectable");
    }

    #[test]
    fn quantum_noise_degrades_recon_and_scales_with_flux() {
        // End-to-end MVCT noise: inject quantum noise into the disk sinogram, then
        // measure the reconstructed interior noise (roi std). Noise must exceed the
        // (near-zero) noiseless recon ripple, and higher photon flux must lower it.
        use crate::noise::add_quantum_noise;
        use helios_analysis::roi_statistics;
        let mu0 = 0.04;
        let phantom = disk_phantom(mu0, 25.0);
        let angles = uniform_angles(180);
        let offsets = uniform_offsets(45.0, 181);
        let clean_sino = parallel_beam_radon(&phantom, &angles, &offsets, 400.0, 0.25);

        let interior = [18usize, 18, 0];
        let interior_hi = [23usize, 23, 1];
        let clean_std = roi_statistics(
            &filtered_back_projection(&clean_sino, &recon_grid()),
            interior,
            interior_hi,
        )
        .std;

        // Low flux → visible noise; high flux → less. Both exceed the clean ripple.
        let low = add_quantum_noise(&clean_sino, 1.0e4, 20260701);
        let high = add_quantum_noise(&clean_sino, 1.0e6, 20260701);
        let std_low = roi_statistics(
            &filtered_back_projection(&low, &recon_grid()),
            interior,
            interior_hi,
        )
        .std;
        let std_high = roi_statistics(
            &filtered_back_projection(&high, &recon_grid()),
            interior,
            interior_hi,
        )
        .std;

        assert!(
            std_low > clean_std,
            "noisy recon std {std_low} not > clean {clean_std}"
        );
        assert!(
            std_low > std_high,
            "lower flux must be noisier: {std_low} !> {std_high}"
        );
    }

    #[test]
    fn ram_lak_kernel_has_expected_structure() {
        let k = ram_lak_kernel::<f64>(4, 0.1); // len 4 → indices n=-3..3
        let base = 3;
        // Even taps (n=±2) are zero.
        assert_eq!(k[base + 2], 0.0);
        assert_eq!(k[base - 2], 0.0);
        // Peak at n=0 is positive; odd taps negative.
        assert!(k[base] > 0.0);
        assert!(k[base + 1] < 0.0 && k[base + 3] < 0.0);
        // Symmetric.
        assert_relative_eq!(k[base + 1], k[base - 1], epsilon = 1e-15);
        assert_relative_eq!(k[base + 3], k[base - 3], epsilon = 1e-15);
    }

    #[test]
    fn reconstruction_is_generic_over_scalar_f32() {
        // f32 disk phantom, radon, FBP — the pipeline is generic over Scalar.
        let n = 161;
        let spacing = 0.5_f32;
        let grid = VoxelGrid::<f32>::axis_aligned(
            [n, n, 1],
            [spacing, spacing, spacing],
            Point3::new(0.0, 0.0, 0.0),
        )
        .unwrap();
        let c = (n - 1) as f32 * spacing / 2.0;
        let phantom = Volume::from_shape_fn(grid, move |idx| {
            let dx = idx[0] as f32 * spacing - c;
            let dy = idx[1] as f32 * spacing - c;
            if (dx * dx + dy * dy).sqrt() <= 25.0 {
                0.04_f32
            } else {
                0.0
            }
        });
        let angles: Vec<f32> = (0..90)
            .map(|a| a as f32 * std::f32::consts::PI / 90.0)
            .collect();
        let offsets: Vec<f32> = (0..121).map(|j| -45.0 + j as f32 * 90.0 / 120.0).collect();
        let sino = parallel_beam_radon(&phantom, &angles, &offsets, 400.0, 0.5);
        let recon = filtered_back_projection(
            &sino,
            &VoxelGrid::<f32>::axis_aligned(
                [41, 41, 1],
                [2.0, 2.0, 2.0],
                Point3::new(0.0, 0.0, 0.0),
            )
            .unwrap(),
        );
        let centre = recon.get(20, 20, 0).unwrap();
        assert_relative_eq!(centre, 0.04_f32, max_relative = 0.2);
    }
}
