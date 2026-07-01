//! SIRT (Simultaneous Iterative Reconstruction Technique).
//!
//! An algebraic MVCT reconstruction: iteratively refine an attenuation estimate
//! so its forward projection matches the measured sinogram. Where filtered
//! back-projection ([`crate::filtered_back_projection`]) is a one-shot analytic
//! inverse, SIRT is iterative and robust to noise and sparse/limited-angle data —
//! the regime where FBP streaks.
//!
//! Update (normalized SIRT): `x ← max(0, x + λ · C⁻¹ ⊙ Aᵀ( R⁻¹ ⊙ (b − A x) ))`,
//! where `A` is the parallel-beam Radon forward projector, `Aᵀ` the back-projector,
//! `R⁻¹ = 1/(A·1)` (per-ray chord length) and `C⁻¹ = 1/(Aᵀ·1)` (per-voxel hit
//! weight). The non-negativity projection encodes `μ ≥ 0`. For consistent data
//! the residual decreases monotonically toward the least-squares solution.

use crate::backproject::back_project_rows;
use crate::radon::{parallel_beam_radon, Sinogram};
use helios_domain::{Volume, VoxelGrid};
use helios_math::{GeometryScalar, NumericElement};

/// Reconstruct the attenuation slice from `sinogram` by `iterations` of SIRT onto
/// the `recon` grid.
///
/// `source_distance_mm`/`step_mm` parameterize the forward projector used each
/// iteration (must match how the sinogram was formed); `relaxation` is the SIRT
/// step `λ` (stable for `0 < λ < 2`; `1.0` is the standard choice). Reconstructed
/// values are non-negative linear attenuation `μ` (cm⁻¹).
#[must_use]
pub fn sirt_reconstruction<T: GeometryScalar>(
    sinogram: &Sinogram<T>,
    recon: &VoxelGrid<T>,
    source_distance_mm: T,
    step_mm: T,
    iterations: usize,
    relaxation: T,
) -> Volume<T> {
    let zero = <T as NumericElement>::ZERO;
    let one = <T as NumericElement>::ONE;
    let angles = sinogram.angles();
    let offsets = sinogram.offsets();
    let (n_ang, n_off) = sinogram.dims();
    let n = n_ang * n_off;

    // Row normalization R⁻¹ = 1/(A·1): forward-project a unit volume → each ray's
    // chord length. Column normalization is the back-projected all-ones sinogram.
    let ones_vol = Volume::from_shape_fn(*recon, |_| one);
    let a_ones = parallel_beam_radon(&ones_vol, angles, offsets, source_distance_mm, step_mm);
    let r_inv: Vec<T> = (0..n)
        .map(|idx| {
            let v = a_ones.get(idx / n_off, idx % n_off);
            if v > zero {
                v.recip()
            } else {
                zero
            }
        })
        .collect();
    let ones_rows = vec![one; n];
    let c_weight = back_project_rows(angles, offsets, &ones_rows, recon, one);

    let mut x = Volume::zeros(*recon);
    for _ in 0..iterations {
        // Row residual R⁻¹ ⊙ (b − A x).
        let ax = parallel_beam_radon(&x, angles, offsets, source_distance_mm, step_mm);
        let mut resid = vec![zero; n];
        for a in 0..n_ang {
            for d in 0..n_off {
                let idx = a * n_off + d;
                resid[idx] = (sinogram.get(a, d) - ax.get(a, d)) * r_inv[idx];
            }
        }
        // Column-normalized, relaxed, non-negativity-projected update.
        let correction = back_project_rows(angles, offsets, &resid, recon, one);
        let x_next = Volume::from_shape_fn(*recon, |idx| {
            let [i, j, k] = idx;
            let c = c_weight.get(i, j, k).expect("in grid");
            let c_inv = if c > zero { c.recip() } else { zero };
            let updated = x.get(i, j, k).expect("in grid")
                + relaxation * c_inv * correction.get(i, j, k).expect("in grid");
            if updated > zero {
                updated
            } else {
                zero
            }
        });
        x = x_next;
    }
    x
}

#[cfg(test)]
mod tests {
    use crate::radon::parallel_beam_radon;
    use crate::sirt_reconstruction;
    use helios_analysis::{roi_statistics, volume_relative_l2_error};
    use helios_domain::{Volume, VoxelGrid};
    use helios_math::Point3;

    fn recon_grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([21, 21, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0)).unwrap()
    }

    // Disk phantom on the reconstruction grid (μ₀ inside radius, 0 outside).
    fn disk(mu0: f64, radius_mm: f64) -> Volume<f64> {
        let grid = recon_grid();
        let centre = 20.0; // (21-1)/2 · 2 mm
        Volume::from_shape_fn(grid, move |idx| {
            let dx = idx[0] as f64 * 2.0 - centre;
            let dy = idx[1] as f64 * 2.0 - centre;
            if (dx * dx + dy * dy).sqrt() <= radius_mm {
                mu0
            } else {
                0.0
            }
        })
    }

    fn angles(n: usize) -> Vec<f64> {
        (0..n)
            .map(|a| a as f64 * std::f64::consts::PI / n as f64)
            .collect()
    }
    fn offsets(half: f64, n: usize) -> Vec<f64> {
        let ds = 2.0 * half / (n - 1) as f64;
        (0..n).map(|j| -half + j as f64 * ds).collect()
    }

    #[test]
    fn sirt_inverts_its_own_forward_model() {
        // Consistent-data (inverse-crime) test of the SOLVER: b = A·phantom, so the
        // phantom is the exact least-squares solution and SIRT must converge to it.
        // (Discretization accuracy vs the analytical sinogram is the FBP tests' job.)
        let phantom = disk(0.04, 14.0);
        let (ang, off) = (angles(40), offsets(24.0, 41));
        let b = parallel_beam_radon(&phantom, &ang, &off, 400.0, 0.5);

        let x = sirt_reconstruction(&b, &recon_grid(), 400.0, 0.5, 30, 1.0);
        // Interior accuracy (the well-determined region): mean μ recovers μ₀ within
        // the same 15% reconstruction tolerance used for FBP — the whole-image L2 is
        // edge-Gibbs-dominated at a sharp disk boundary, so it is a looser sanity
        // bound here (both FBP and SIRT show this at the discontinuity).
        let interior = roi_statistics(&x, [8, 8, 0], [13, 13, 1]);
        assert!(
            (interior.mean - 0.04).abs() / 0.04 < 0.15,
            "interior mean {} not within 15% of μ₀",
            interior.mean
        );
        let err = volume_relative_l2_error(&x, &phantom).unwrap();
        assert!(err < 0.2, "SIRT relative L2 error {err} should be < 20%");
    }

    #[test]
    fn more_iterations_reduce_the_error() {
        let phantom = disk(0.04, 14.0);
        let (ang, off) = (angles(40), offsets(24.0, 41));
        let b = parallel_beam_radon(&phantom, &ang, &off, 400.0, 0.5);
        let grid = recon_grid();
        let err5 = volume_relative_l2_error(
            &sirt_reconstruction(&b, &grid, 400.0, 0.5, 5, 1.0),
            &phantom,
        )
        .unwrap();
        let err25 = volume_relative_l2_error(
            &sirt_reconstruction(&b, &grid, 400.0, 0.5, 25, 1.0),
            &phantom,
        )
        .unwrap();
        assert!(
            err25 < err5,
            "error must fall with iterations: {err25} !< {err5}"
        );
    }

    #[test]
    fn zero_sinogram_reconstructs_zero() {
        let (ang, off) = (angles(20), offsets(24.0, 21));
        let empty = parallel_beam_radon(&Volume::zeros(recon_grid()), &ang, &off, 400.0, 1.0);
        let x = sirt_reconstruction(&empty, &recon_grid(), 400.0, 1.0, 10, 1.0);
        assert!((0..21).all(|i| (0..21).all(|j| x.get(i, j, 0).unwrap().abs() < 1e-12)));
    }

    #[test]
    fn sirt_is_generic_over_scalar_f32() {
        let grid = VoxelGrid::<f32>::axis_aligned(
            [15, 15, 1],
            [2.0, 2.0, 2.0],
            Point3::new(0.0, 0.0, 0.0),
        )
        .unwrap();
        let phantom = Volume::from_shape_fn(grid, |idx| {
            let (dx, dy) = (idx[0] as f32 * 2.0 - 14.0, idx[1] as f32 * 2.0 - 14.0);
            if (dx * dx + dy * dy).sqrt() <= 10.0 {
                0.04
            } else {
                0.0
            }
        });
        let ang: Vec<f32> = (0..30)
            .map(|a| a as f32 * std::f32::consts::PI / 30.0)
            .collect();
        let off: Vec<f32> = (0..31).map(|j| -20.0 + j as f32 * 40.0 / 30.0).collect();
        let b = parallel_beam_radon(&phantom, &ang, &off, 400.0, 0.5);
        let x = sirt_reconstruction(&b, &grid, 400.0, 0.5, 30, 1.0);
        // Interior mean recovers μ₀ within the 15% reconstruction tolerance.
        let interior = roi_statistics(&x, [5, 5, 0], [10, 10, 1]);
        assert!(
            (interior.mean - 0.04_f32).abs() / 0.04 < 0.15,
            "f32 interior mean {} not within 15% of μ₀",
            interior.mean
        );
    }
}
