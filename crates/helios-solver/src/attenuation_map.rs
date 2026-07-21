//! CT-number → linear-attenuation map (deterministic material-property engine).

use helios_domain::Volume;
use helios_math::Scalar;
use helios_physics::{mass_density_from_hu, MassAttenuation};

/// Map a CT volume (Hounsfield units) to a linear-attenuation volume `μ` (cm⁻¹)
/// at a fixed photon energy.
///
/// For each voxel: `ρ = mass_density_from_hu(HU, ρ_water)` and
/// `μ = (μ/ρ)·ρ`, using the water mass-attenuation coefficient `mass_attenuation`
/// at the beam energy. This is the **Compton-dominated MV approximation** (the
/// TomoTherapy 6 MV regime): at MeV energies photon attenuation scales with
/// electron ≈ mass density at an approximately material-independent `(μ/ρ)`, so a
/// single water `(μ/ρ)` scaled by voxel density is an accurate first-order model.
/// (A kV/energy-dependent, material-segmented model is a later refinement.)
///
/// The output volume shares the input's [`VoxelGrid`](helios_domain::VoxelGrid).
/// All voxel values are non-negative (density is clamped at zero below air).
#[must_use]
pub fn attenuation_map<T: Scalar>(
    ct_hu: &Volume<T>,
    mass_attenuation: MassAttenuation<T>,
    water_density_g_cm3: T,
) -> Volume<T> {
    let grid = *ct_hu.grid();
    let mu_over_rho = mass_attenuation.get();
    Volume::from_shape_fn(grid, |idx| {
        let hu = ct_hu
            .get(idx[0], idx[1], idx[2])
            .expect("from_shape_fn iterates indices within the shared grid");
        let density = mass_density_from_hu(hu, water_density_g_cm3);
        mu_over_rho * density
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    fn grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([3, 4, 5], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .expect("valid grid")
    }

    fn water_mass_attenuation() -> MassAttenuation<f64> {
        // Representative water μ/ρ magnitude (cm²/g); the engine is verified by
        // the defining relation, not by this specific value.
        MassAttenuation::new(0.06).expect("valid coefficient")
    }

    #[test]
    fn uniform_water_maps_to_constant_mu() {
        // All HU = 0 (water) → μ = (μ/ρ)·ρ_water = 0.06·1.0 everywhere.
        let ct = Volume::from_shape_fn(grid(), |_| 0.0);
        let mu = attenuation_map(&ct, water_mass_attenuation(), 1.0);
        for i in 0..3 {
            for j in 0..4 {
                for k in 0..5 {
                    assert_relative_eq!(mu.get(i, j, k).unwrap(), 0.06, epsilon = 1e-15);
                }
            }
        }
    }

    #[test]
    fn air_maps_to_zero_and_bone_scales_with_density() {
        // Air (−1000 HU) → ρ=0 → μ=0; HU=1000 → ρ=2 → μ=0.12.
        let ct = Volume::from_shape_vec(
            VoxelGrid::axis_aligned([2, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap(),
            vec![-1000.0, 1000.0],
        )
        .unwrap();
        let mu = attenuation_map(&ct, water_mass_attenuation(), 1.0);
        assert_relative_eq!(mu.get(0, 0, 0).unwrap(), 0.0, epsilon = 1e-15);
        assert_relative_eq!(mu.get(1, 0, 0).unwrap(), 0.12, epsilon = 1e-15);
    }

    #[test]
    fn engine_matches_direct_per_voxel_formula() {
        // Differential oracle over a heterogeneous HU field: engine == closed form.
        let ct = Volume::from_shape_fn(grid(), |idx| {
            // A varied HU field spanning air→soft-tissue→bone.
            -800.0 + 120.0 * idx[0] as f64 + 90.0 * idx[1] as f64 + 60.0 * idx[2] as f64
        });
        let mass_atten = water_mass_attenuation();
        let water_density = 1.0;
        let mu = attenuation_map(&ct, mass_atten, water_density);
        for i in 0..3 {
            for j in 0..4 {
                for k in 0..5 {
                    let hu = ct.get(i, j, k).unwrap();
                    let expected = mass_atten.get() * mass_density_from_hu(hu, water_density);
                    assert_relative_eq!(mu.get(i, j, k).unwrap(), expected, epsilon = 1e-15);
                }
            }
        }
    }

    #[test]
    fn output_grid_matches_input() {
        let ct = Volume::from_shape_fn(grid(), |_| 0.0);
        let mu = attenuation_map(&ct, water_mass_attenuation(), 1.0);
        assert_eq!(mu.grid().dims(), ct.grid().dims());
    }

    #[test]
    fn engine_is_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([2, 2, 2], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let ct = Volume::from_shape_fn(g, |_| 0.0_f32);
        let mu = attenuation_map(&ct, MassAttenuation::new(0.06_f32).unwrap(), 1.0);
        assert_relative_eq!(mu.get(0, 0, 0).unwrap(), 0.06_f32, epsilon = 1e-6);
    }
}
