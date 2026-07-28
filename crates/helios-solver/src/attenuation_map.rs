//! CT-number → linear-attenuation map (deterministic material-property engine).

use core::fmt;

use aequitas::systems::si::{
    quantities::MassDensity as DensityQuantity,
    units::{GramPerCubicCentimeter, PerCentimeter},
};
use eunomia::UnitScalar;
use helios_domain::Volume;
use helios_math::Scalar;
use helios_physics::mass_density_from_hu;
use hyperion::{TransportError, coefficient::MassAttenuation};
use proteus::{InvalidProperty, MassDensity};

/// Failure while converting a CT volume into linear attenuation.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum AttenuationMapError<T> {
    /// A calibrated voxel density violated the material-property contract.
    InvalidDensity(InvalidProperty<T>),
    /// Mass-to-linear attenuation evaluation failed.
    Transport(TransportError<T>),
}

impl<T> From<InvalidProperty<T>> for AttenuationMapError<T> {
    fn from(error: InvalidProperty<T>) -> Self {
        Self::InvalidDensity(error)
    }
}

impl<T> From<TransportError<T>> for AttenuationMapError<T> {
    fn from(error: TransportError<T>) -> Self {
        Self::Transport(error)
    }
}

impl<T: fmt::Debug> fmt::Display for AttenuationMapError<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDensity(error) => write!(formatter, "invalid calibrated density: {error}"),
            Self::Transport(error) => write!(formatter, "attenuation conversion failed: {error}"),
        }
    }
}

impl<T: fmt::Debug> core::error::Error for AttenuationMapError<T> {}

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
///
/// # Errors
///
/// Returns [`AttenuationMapError::InvalidDensity`] if CT calibration produces a
/// non-finite or negative density, and [`AttenuationMapError::Transport`] if the
/// mass-to-linear attenuation product is non-finite.
pub fn attenuation_map<T: Scalar + UnitScalar>(
    ct_hu: &Volume<T>,
    mass_attenuation: MassAttenuation<T>,
    water_density: DensityQuantity<T>,
) -> Result<Volume<T>, AttenuationMapError<T>> {
    let water_density_g_cm3: T = water_density.in_unit::<GramPerCubicCentimeter>();
    let grid = *ct_hu.grid();
    let values = ct_hu
        .as_slice()
        .iter()
        .copied()
        .map(|hu| {
            let density = mass_density_from_hu(hu, water_density_g_cm3);
            let density = MassDensity::new(DensityQuantity::from_unit::<GramPerCubicCentimeter>(
                density,
            ))?;
            Ok(mass_attenuation
                .to_linear(density)?
                .in_unit::<PerCentimeter>())
        })
        .collect::<Result<Vec<_>, AttenuationMapError<T>>>()?;
    Ok(Volume::from_shape_vec(grid, values)
        .expect("invariant: one attenuation value is produced for every CT voxel"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::{quantities::AreaPerMass, units::SquareCentimeterPerGram};
    use eunomia::assert_relative_eq;
    use helios_math::ShippedScalar;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    fn grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([3, 4, 5], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .expect("valid grid")
    }

    fn water_mass_attenuation() -> MassAttenuation<f64> {
        // Representative water μ/ρ magnitude (cm²/g); the engine is verified by
        // the defining relation, not by this specific value.
        MassAttenuation::new(AreaPerMass::from_unit::<SquareCentimeterPerGram>(0.06))
            .expect("valid coefficient")
    }

    #[test]
    fn uniform_water_maps_to_constant_mu() {
        // All HU = 0 (water) → μ = (μ/ρ)·ρ_water = 0.06·1.0 everywhere.
        let ct = Volume::from_shape_fn(grid(), |_| 0.0);
        let mu = attenuation_map(
            &ct,
            water_mass_attenuation(),
            DensityQuantity::from_unit::<GramPerCubicCentimeter>(1.0),
        )
        .expect("fixture calibration is finite");
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
        let mu = attenuation_map(
            &ct,
            water_mass_attenuation(),
            DensityQuantity::from_unit::<GramPerCubicCentimeter>(1.0),
        )
        .expect("fixture calibration is finite");
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
        let water_density = DensityQuantity::from_unit::<GramPerCubicCentimeter>(1.0);
        let mu =
            attenuation_map(&ct, mass_atten, water_density).expect("fixture calibration is finite");
        for i in 0..3 {
            for j in 0..4 {
                for k in 0..5 {
                    let hu = ct.get(i, j, k).unwrap();
                    let expected = mass_atten.in_unit::<SquareCentimeterPerGram>()
                        * mass_density_from_hu(
                            hu,
                            water_density.in_unit::<GramPerCubicCentimeter>(),
                        );
                    assert_relative_eq!(mu.get(i, j, k).unwrap(), expected, epsilon = 1e-15);
                }
            }
        }
    }

    #[test]
    fn output_grid_matches_input() {
        let ct = Volume::from_shape_fn(grid(), |_| 0.0);
        let mu = attenuation_map(
            &ct,
            water_mass_attenuation(),
            DensityQuantity::from_unit::<GramPerCubicCentimeter>(1.0),
        )
        .expect("fixture calibration is finite");
        assert_eq!(mu.grid().dims(), ct.grid().dims());
    }

    /// `mu = (mu/rho) * rho` at unit density is one product plus the HU-to-density
    /// calibration of a zero CT value, so at most a few roundings separate the
    /// result from the mass-attenuation input. Four ulps of `T` bounds that.
    const MASS_TO_LINEAR_ULPS: f64 = 4.0;

    /// Asserts mass-to-linear attenuation conversion in one scalar width.
    fn attenuation_map_scales_mass_coefficient_by_density<T: ShippedScalar>() {
        let zero = T::from_f64(0.0);
        let unit = T::from_f64(1.0);
        let grid = VoxelGrid::<T>::axis_aligned(
            [2, 2, 2],
            [unit; 3],
            Point3::new(zero, zero, zero),
        )
        .expect("valid axis-aligned grid");

        let mass_coefficient = T::from_f64(0.06);
        let ct = Volume::from_shape_fn(grid, |_| zero);
        let coefficient = MassAttenuation::new(AreaPerMass::from_unit::<
            SquareCentimeterPerGram,
        >(mass_coefficient))
        .expect("positive mass-attenuation coefficient");
        let mu = attenuation_map(
            &ct,
            coefficient,
            DensityQuantity::<T>::from_unit::<GramPerCubicCentimeter>(unit),
        )
        .expect("fixture calibration is finite");

        // At 1 g/cm^3 the linear coefficient equals the mass coefficient.
        assert_relative_eq!(
            mu.get(0, 0, 0).unwrap(),
            mass_coefficient,
            max_relative = T::EPSILON * T::from_f64(MASS_TO_LINEAR_ULPS)
        );
    }

    #[test]
    fn attenuation_map_scales_mass_coefficient_by_density_in_single_precision() {
        attenuation_map_scales_mass_coefficient_by_density::<f32>();
    }

    #[test]
    fn attenuation_map_scales_mass_coefficient_by_density_in_double_precision() {
        attenuation_map_scales_mass_coefficient_by_density::<f64>();
    }

    #[test]
    fn invalid_reference_density_preserves_the_proteus_error() {
        let ct = Volume::from_shape_fn(grid(), |_| 0.0);
        let error = attenuation_map(
            &ct,
            water_mass_attenuation(),
            DensityQuantity::from_unit::<GramPerCubicCentimeter>(-1.0_f64),
        )
        .expect_err("negative calibrated density must be rejected");
        assert!(matches!(error, AttenuationMapError::InvalidDensity(_)));
    }
}
