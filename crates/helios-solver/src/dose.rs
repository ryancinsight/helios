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

use aequitas::systems::si::{
    quantities::{Length, ReciprocalLength},
    units::{Centimeter, PerCentimeter},
};
use eunomia::UnitScalar;
use helios_core::constants::{CM_PER_M, MM_PER_CM};
use helios_domain::{Volume, VoxelGrid};
use helios_math::Scalar;
use hyperion::{
    coefficient::{InteractionCoefficient, LinearAttenuation},
    quantity::{OpticalDepth, PathLength},
    TransportError,
};

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
///
/// # Errors
///
/// Returns [`TransportError`] when an attenuation coefficient is negative or
/// non-finite, or when an optical-depth product or partial sum is non-finite.
pub fn primary_fluence_parallel_x<T: Scalar + UnitScalar>(
    mu: &Volume<T>,
    incident_fluence: T,
) -> Result<Volume<T>, TransportError<T>> {
    let grid: VoxelGrid<T> = *mu.grid();
    // Voxel spacing along the beam axis, in cm.
    let sx_cm = grid.spacing()[0] * T::from_f64(MM_PER_CM).recip();
    let path = PathLength::new(Length::from_unit::<Centimeter>(sx_cm))?;
    let [nx, ny, nz] = grid.dims();
    let mut values = Vec::with_capacity(grid.num_voxels());
    let mut optical_depths = vec![OpticalDepth::zero(); ny * nz];

    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                // Optical depth from the entry face to the centre of column i:
                // Σ_{i'<i} μ(i',j,k) · sx. Hyperion owns coefficient validation,
                // dimensional multiplication, finite reduction, and transmission.
                let ray_index = j * nz + k;
                let coefficient = InteractionCoefficient::<T, LinearAttenuation>::new(
                    ReciprocalLength::from_unit::<PerCentimeter>(
                        mu.get(i, j, k).expect("column index within grid"),
                    ),
                )?;
                values.push(
                    incident_fluence
                        * optical_depths[ray_index]
                            .transmission()
                            .into_quantity()
                            .into_base(),
                );
                if i + 1 < nx {
                    optical_depths[ray_index] =
                        optical_depths[ray_index].checked_add(coefficient.optical_depth(path)?)?;
                }
            }
        }
    }

    Ok(Volume::from_shape_vec(grid, values)
        .expect("invariant: one primary-fluence value is produced for every voxel"))
}

/// Build a normalized forward (downstream) exponential dose-deposition kernel.
///
/// `kernel[d] ∝ exp(-d·voxel_spacing / range)` for `d ∈ [0, taps)`, normalized so
/// `Σ kernel = 1` (energy-conserving). `range` is the characteristic deposition
/// range (electron transport scale); `voxel_spacing` is the voxel spacing along
/// the beam. Both values retain their Aequitas units through this API and are
/// converted to centimetres only at the exponential formula boundary. Returns
/// an empty vector when `taps == 0`.
#[must_use]
pub fn exponential_deposition_kernel<T: Scalar>(
    range: Length<T>,
    voxel_spacing: Length<T>,
    taps: usize,
) -> Vec<T> {
    let mut kernel = Vec::with_capacity(taps);
    let mut sum = <T as helios_math::NumericElement>::ZERO;
    let range_cm = range.into_base() * T::from_f64(CM_PER_M);
    let voxel_cm = voxel_spacing.into_base() * T::from_f64(CM_PER_M);
    let inv_range = range_cm.recip();
    for d in 0..taps {
        let distance = T::from_f64(d as f64) * voxel_cm;
        let weight = (-(distance * inv_range)).exp();
        kernel.push(weight);
        sum += weight;
    }
    if sum > <T as helios_math::NumericElement>::ZERO {
        let inv_sum = sum.recip();
        for w in &mut kernel {
            *w *= inv_sum;
        }
    }
    kernel
}

/// Dose by 1-D convolution-superposition of a TERMA/energy-release volume with a
/// forward dose-deposition `kernel` along the beam axis (**+x**).
///
/// Energy released at column `i'` is deposited downstream at columns `i ≥ i'`:
///
/// ```text
/// D(i) = Σ_{d=0}^{min(i, taps-1)} TERMA(i-d) · kernel[d]
/// ```
///
/// This is the discrete collapsed-cone / convolution-superposition step: with a
/// spread-out kernel it reproduces the photon depth-dose **build-up** (shallow
/// voxels receive less than the downstream maximum). A single-tap `[1]` kernel is
/// the identity (dose = TERMA); a normalized kernel conserves energy in the
/// interior. Returns zeros for an empty kernel.
#[must_use]
pub fn dose_convolution_x<T: Scalar>(terma: &Volume<T>, kernel: &[T]) -> Volume<T> {
    let grid: VoxelGrid<T> = *terma.grid();
    let taps = kernel.len();
    Volume::from_shape_fn(grid, |idx| {
        let [i, j, k] = idx;
        if taps == 0 {
            return <T as helios_math::NumericElement>::ZERO;
        }
        let max_d = i.min(taps - 1);
        let mut dose = <T as helios_math::NumericElement>::ZERO;
        for (d, &weight) in kernel.iter().enumerate().take(max_d + 1) {
            dose += terma.get(i - d, j, k).expect("column index within grid") * weight;
        }
        dose
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_math::Point3;
    use helios_math::ShippedScalar;

    fn grid() -> VoxelGrid<f64> {
        // 2 mm spacing along x → 0.2 cm per column.
        VoxelGrid::axis_aligned([6, 2, 2], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid")
    }

    fn length_cm(value: f64) -> Length<f64> {
        Length::from_base(value * 0.01)
    }

    #[test]
    fn homogeneous_medium_gives_exponential_depth_curve() {
        // Uniform μ = 0.3 cm⁻¹, Ψ₀ = 5.0. Ψ(i) = 5·exp(-0.3 · i·0.2).
        let mu = Volume::from_shape_fn(grid(), |_| 0.3);
        let psi = primary_fluence_parallel_x(&mu, 5.0).expect("valid attenuation volume");
        for i in 0..6 {
            let depth_cm = i as f64 * 0.2;
            let expected = 5.0 * (-0.3 * depth_cm).exp();
            assert_relative_eq!(psi.get(i, 0, 0).unwrap(), expected, epsilon = 1e-12);
        }
    }

    #[test]
    fn entry_column_is_unattenuated() {
        let mu = Volume::from_shape_fn(grid(), |_| 0.9);
        let psi = primary_fluence_parallel_x(&mu, 2.0).expect("valid attenuation volume");
        assert_relative_eq!(psi.get(0, 1, 1).unwrap(), 2.0, epsilon = 1e-15);
    }

    #[test]
    fn heterogeneous_columns_accumulate_optical_depth() {
        // μ varies along x: column i has μ = 0.1·(i+1). τ to column i =
        // Σ_{i'<i} 0.1·(i'+1) · 0.2.
        let mu = Volume::from_shape_fn(grid(), |idx| 0.1 * (idx[0] as f64 + 1.0));
        let psi = primary_fluence_parallel_x(&mu, 1.0).expect("valid attenuation volume");
        let mut tau = 0.0_f64;
        for i in 0..6 {
            let expected = (-tau).exp();
            assert_relative_eq!(psi.get(i, 0, 0).unwrap(), expected, epsilon = 1e-12);
            tau += 0.1 * (i as f64 + 1.0) * 0.2;
        }
    }

    #[test]
    fn delta_kernel_is_identity() {
        // A single unit tap deposits all energy locally → dose == TERMA.
        let terma = Volume::from_shape_fn(grid(), |idx| 1.0 + idx[0] as f64);
        let dose = dose_convolution_x(&terma, &[1.0]);
        for i in 0..6 {
            assert_relative_eq!(
                dose.get(i, 0, 0).unwrap(),
                terma.get(i, 0, 0).unwrap(),
                epsilon = 1e-15
            );
        }
    }

    #[test]
    fn normalized_kernel_conserves_uniform_terma_in_interior() {
        // Uniform TERMA, normalized 3-tap kernel: once all taps fit (i ≥ 2), dose
        // returns to the uniform value (energy conservation); shallow voxels get
        // less (build-up).
        let terma = Volume::from_shape_fn(grid(), |_| 4.0);
        let kernel = [0.5, 0.3, 0.2]; // sums to 1
        let dose = dose_convolution_x(&terma, &kernel);
        assert_relative_eq!(dose.get(0, 0, 0).unwrap(), 4.0 * 0.5, epsilon = 1e-15);
        assert_relative_eq!(dose.get(1, 0, 0).unwrap(), 4.0 * 0.8, epsilon = 1e-15);
        for i in 2..6 {
            assert_relative_eq!(dose.get(i, 0, 0).unwrap(), 4.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn exponential_kernel_is_normalized() {
        let kernel = exponential_deposition_kernel(length_cm(0.5), length_cm(0.2), 8);
        let sum: f64 = kernel.iter().sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-15);
        assert_eq!(kernel.len(), 8);
    }

    #[test]
    fn convolution_produces_physical_buildup() {
        // Exponential TERMA through a spread kernel exhibits photon build-up:
        // the surface voxel receives less than a downstream maximum, then dose
        // falls off. This is the clinical-realism depth-dose signature.
        let g = VoxelGrid::axis_aligned([40, 1, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let mu = Volume::from_shape_fn(g, |_| 0.3);
        let terma = primary_fluence_parallel_x(&mu, 1.0).expect("valid attenuation volume");
        let kernel = exponential_deposition_kernel(length_cm(1.0), length_cm(0.2), 20);
        let dose = dose_convolution_x(&terma, &kernel);

        // Find the peak; it must be deeper than the surface (build-up) and the
        // dose must fall off beyond it.
        let vals: Vec<f64> = (0..40).map(|i| dose.get(i, 0, 0).unwrap()).collect();
        let (peak_i, &peak) = vals
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap();
        assert!(
            peak_i > 0,
            "peak must be below the surface (build-up region)"
        );
        assert!(vals[0] < peak, "surface dose must be below the peak");
        assert!(vals[39] < peak, "dose must fall off beyond the peak");
    }

    #[test]
    fn empty_kernel_gives_zero_dose() {
        let terma = Volume::from_shape_fn(grid(), |_| 3.0);
        let dose = dose_convolution_x(&terma, &[]);
        for i in 0..6 {
            assert_relative_eq!(dose.get(i, 0, 0).unwrap(), 0.0, epsilon = 1e-15);
        }
    }

    /// Relative tolerance, in ulps of `T`, for a Beer-Lambert transmission.
    ///
    /// The value accumulates optical depth across three voxels, converts spacing
    /// to centimetres once, evaluates `exp`, then scales by the incident
    /// fluence. Each arithmetic step contributes at most one rounding of
    /// relative size `T::EPSILON`, and `exp` is not required by IEEE 754 to be
    /// correctly rounded, so it is budgeted at several ulps rather than one.
    /// Sixteen bounds that chain on both shipped widths — 1.9e-6 for `f32`,
    /// 3.6e-15 for `f64` — and is tighter than the 1e-5 absolute bound the
    /// single-precision-only predecessor used, so it cannot mask a divergence
    /// the old assertion would have caught.
    const TRANSMISSION_ULPS: f64 = 16.0;

    /// Asserts the analytical Beer-Lambert transmission in one scalar width.
    ///
    /// Written once and instantiated per shipped scalar type below: a generic
    /// kernel verified at a single concrete type leaves every other
    /// monomorphization callers can build unverified.
    fn primary_fluence_matches_beer_lambert<T: ShippedScalar>() {
        let zero = T::from_f64(0.0);
        let spacing_mm = T::from_f64(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([4, 1, 1], [spacing_mm; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        let mu = T::from_f64(0.3);
        let incident = T::from_f64(5.0);
        let volume = Volume::from_shape_fn(grid, |_| mu);
        let psi = primary_fluence_parallel_x(&volume, incident).expect("valid attenuation volume");

        // Column 3 sits behind three voxels: 3 x 2 mm = 0.6 cm of path.
        let depth_cm = T::from_f64(3.0) * spacing_mm * T::from_f64(MM_PER_CM).recip();
        let expected = incident * (-(mu * depth_cm)).exp();

        assert_relative_eq!(
            psi.get(3, 0, 0).unwrap(),
            expected,
            max_relative = T::EPSILON * T::from_f64(TRANSMISSION_ULPS)
        );
    }

    #[test]
    fn primary_fluence_matches_beer_lambert_in_single_precision() {
        primary_fluence_matches_beer_lambert::<f32>();
    }

    #[test]
    fn primary_fluence_matches_beer_lambert_in_double_precision() {
        primary_fluence_matches_beer_lambert::<f64>();
    }

    #[test]
    fn negative_attenuation_is_rejected_before_transmission() {
        let mu = Volume::from_shape_fn(grid(), |[i, _, _]| if i == 5 { -0.1 } else { 0.3 });
        let error = primary_fluence_parallel_x(&mu, 1.0).expect_err("negative mu must fail");
        assert!(matches!(
            error,
            TransportError::InvalidValue {
                field: hyperion::ValueKind::LinearAttenuation,
                ..
            }
        ));
    }
}
