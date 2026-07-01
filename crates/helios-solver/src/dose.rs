//! Primary photon transport — the first stage of the deterministic dose engine.
//!
//! A collapsed-cone / convolution-superposition dose calculation proceeds in two
//! stages: (1) **primary transport** — attenuate the incident energy fluence
//! through the patient to get the primary fluence (and TERMA) at every voxel; then
//! (2) **kernel superposition** — convolve TERMA with a dose-deposition kernel to
//! spread energy to dose. This module implements stage (1) for a parallel beam;
//! stage (2) (scatter kernel) is tracked separately.
//!
//! For a parallel beam entering along +x, the primary energy fluence at depth is
//! the Beer–Lambert–attenuated incident fluence:
//!
//! ```text
//! Ψ(x) = Ψ₀ · exp(−∫₀ˣ μ dl)
//! ```
//!
//! In a homogeneous medium this is the exact exponential `Ψ₀·exp(−μx)` — the
//! analytical oracle used in the tests.

use helios_domain::{Volume, VoxelGrid};
use helios_math::Scalar;

/// Millimetres per centimetre — grid spacing is mm, `μ` is cm⁻¹.
const MM_PER_CM: f64 = 10.0;

/// Attenuated primary energy fluence for a parallel beam entering along **+x**.
///
/// `mu` is the linear-attenuation volume (cm⁻¹); `incident_fluence` is `Ψ₀` at the
/// entry face (x = voxel column 0). Returns a volume of primary energy fluence
/// `Ψ` at each voxel centre, accumulating optical depth column-by-column so the
/// cost is one pass over the grid.
///
/// The centre of column `i` is treated as one voxel-spacing (`sx`) further into
/// the medium than column `i−1` (centre-to-centre attenuation), so column 0 sees
/// the unattenuated `Ψ₀`.
#[must_use]
pub fn primary_fluence_parallel_x<T: Scalar>(mu: &Volume<T>, incident_fluence: T) -> Volume<T> {
    let grid: VoxelGrid<T> = *mu.grid();
    // Voxel spacing along the beam axis, in cm.
    let sx_cm = grid.spacing()[0] * T::from_f64(MM_PER_CM).recip();

    Volume::from_shape_fn(grid, |idx| {
        let [i, j, k] = idx;
        // Optical depth from the entry face to the centre of column i:
        // Σ_{i'<i} μ(i',j,k) · sx.
        let mut tau = <T as helios_math::NumericElement>::ZERO;
        for col in 0..i {
            tau += mu.get(col, j, k).expect("column index within grid") * sx_cm;
        }
        incident_fluence * (-tau).exp()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use helios_math::Point3;

    fn grid() -> VoxelGrid<f64> {
        // 2 mm spacing along x → 0.2 cm per column.
        VoxelGrid::axis_aligned([6, 2, 2], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid")
    }

    #[test]
    fn homogeneous_medium_gives_exponential_depth_curve() {
        // Uniform μ = 0.3 cm⁻¹, Ψ₀ = 5.0. Ψ(i) = 5·exp(-0.3 · i·0.2).
        let mu = Volume::from_shape_fn(grid(), |_| 0.3);
        let psi = primary_fluence_parallel_x(&mu, 5.0);
        for i in 0..6 {
            let depth_cm = i as f64 * 0.2;
            let expected = 5.0 * (-0.3 * depth_cm).exp();
            assert_relative_eq!(psi.get(i, 0, 0).unwrap(), expected, epsilon = 1e-12);
        }
    }

    #[test]
    fn entry_column_is_unattenuated() {
        let mu = Volume::from_shape_fn(grid(), |_| 0.9);
        let psi = primary_fluence_parallel_x(&mu, 2.0);
        assert_relative_eq!(psi.get(0, 1, 1).unwrap(), 2.0, epsilon = 1e-15);
    }

    #[test]
    fn heterogeneous_columns_accumulate_optical_depth() {
        // μ varies along x: column i has μ = 0.1·(i+1). τ to column i =
        // Σ_{i'<i} 0.1·(i'+1) · 0.2.
        let mu = Volume::from_shape_fn(grid(), |idx| 0.1 * (idx[0] as f64 + 1.0));
        let psi = primary_fluence_parallel_x(&mu, 1.0);
        let mut tau = 0.0_f64;
        for i in 0..6 {
            let expected = (-tau).exp();
            assert_relative_eq!(psi.get(i, 0, 0).unwrap(), expected, epsilon = 1e-12);
            tau += 0.1 * (i as f64 + 1.0) * 0.2;
        }
    }

    #[test]
    fn primary_fluence_is_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([4, 1, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(g, |_| 0.3_f32);
        let psi = primary_fluence_parallel_x(&mu, 5.0_f32);
        let expected = 5.0_f32 * (-0.3_f32 * 3.0 * 0.2).exp(); // column 3, depth 0.6 cm
        assert_relative_eq!(psi.get(3, 0, 0).unwrap(), expected, epsilon = 1e-5);
    }
}
