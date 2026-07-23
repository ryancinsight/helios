//! Gamma index: combined dose-difference / distance-to-agreement comparison.
//!
//! Low's gamma (Med. Phys. 25, 1998) compares an evaluated dose distribution to a
//! reference. For each reference point `r`,
//!
//! ```text
//! γ(r) = min_e sqrt( |x_e − x_r|² / Δd²  +  (D_e − D_r)² / ΔD² )
//! ```
//!
//! where `Δd` is the distance-to-agreement (mm) and `ΔD` the dose-difference
//! criterion. `γ ≤ 1` passes. With a 3%/2 mm criterion, `Δd = 2 mm` and
//! `ΔD = 0.03 · D_norm` (global normalization).

use aequitas::systems::si::quantities::{AbsorbedDose, Length};
use helios_core::HeliosError;
use helios_domain::Volume;
use helios_math::{NumericElement, Scalar};

fn require_positive_finite<T: Scalar>(value: T, field: &'static str) -> Result<(), HeliosError> {
    if !value.is_finite() || value <= <T as NumericElement>::ZERO {
        return Err(HeliosError::InvalidDomainValue {
            field,
            value: value.to_f64(),
            reason: "must be finite and strictly positive",
        });
    }
    Ok(())
}

/// Compute the 3-D gamma index of `evaluated` against `reference` on a shared grid.
///
/// - `dose_diff_criterion`: fractional dose-difference criterion (e.g. `0.03`).
/// - `dta`: distance-to-agreement (e.g. `Length::from_base(2.0)`).
/// - `normalization_dose`: global normalization dose for `ΔD` (e.g. the reference
///   maximum or the prescription).
/// - `search_radius`: neighborhood radius searched for the minimizing point;
///   choose ≥ a few × `dta` so the true minimum is not truncated.
///
/// Returns a gamma value per reference voxel.
///
/// # Errors
/// Returns [`HeliosError`] if the two volumes do not share an identical grid, or
/// if any criterion is non-finite or non-positive.
pub fn gamma_index_3d<T: Scalar>(
    reference: &Volume<T>,
    evaluated: &Volume<T>,
    dose_diff_criterion: T,
    dta: Length<T>,
    normalization_dose: AbsorbedDose<T>,
    search_radius: Length<T>,
) -> Result<Volume<T>, HeliosError> {
    gamma_impl(
        reference,
        evaluated,
        dose_diff_criterion,
        dta,
        Norm::Global(normalization_dose),
        search_radius,
    )
}

/// **Local-normalization** gamma index: the dose-difference criterion `ΔD` is a
/// fraction of the **local reference dose** at each point (`ΔD(r) = criterion·D_r`)
/// rather than a single global value.
///
/// Local normalization is stricter in low-dose regions (where global normalization
/// is lenient), the appropriate choice when relative agreement everywhere matters.
/// Reference points below `low_dose_cutoff` are excluded (their gamma is set to `0`)
/// — the standard low-dose threshold that also avoids dividing by a vanishing `ΔD`.
///
/// # Errors
/// As [`gamma_index_3d`], with `low_dose_cutoff` required finite and positive.
pub fn gamma_index_3d_local<T: Scalar>(
    reference: &Volume<T>,
    evaluated: &Volume<T>,
    dose_diff_criterion: T,
    dta: Length<T>,
    low_dose_cutoff: AbsorbedDose<T>,
    search_radius: Length<T>,
) -> Result<Volume<T>, HeliosError> {
    gamma_impl(
        reference,
        evaluated,
        dose_diff_criterion,
        dta,
        Norm::Local {
            cutoff: low_dose_cutoff,
        },
        search_radius,
    )
}

/// Dose-difference normalization mode for the gamma index.
#[derive(Debug, Clone, Copy)]
enum Norm<T: Scalar> {
    /// Global: `ΔD = criterion · dose` for a single normalization `dose`.
    Global(AbsorbedDose<T>),
    /// Local: `ΔD = criterion · D_r`; points below `cutoff` are excluded.
    Local { cutoff: AbsorbedDose<T> },
}

/// Shared gamma computation for [`gamma_index_3d`] (global) and
/// [`gamma_index_3d_local`] (local), selected by `norm`.
fn gamma_impl<T: Scalar>(
    reference: &Volume<T>,
    evaluated: &Volume<T>,
    dose_diff_criterion: T,
    dta: Length<T>,
    norm: Norm<T>,
    search_radius: Length<T>,
) -> Result<Volume<T>, HeliosError> {
    if reference.grid() != evaluated.grid() {
        return Err(HeliosError::InvalidDomainValue {
            field: "gamma_index_3d::grid",
            value: f64::NAN,
            reason: "reference and evaluated volumes must share an identical grid",
        });
    }
    require_positive_finite(dose_diff_criterion, "gamma::dose_diff_criterion")?;
    require_positive_finite(*dta.as_base(), "gamma::dta")?;
    require_positive_finite(*search_radius.as_base(), "gamma::search_radius")?;
    match norm {
        Norm::Global(d) => require_positive_finite(*d.as_base(), "gamma::normalization_dose")?,
        Norm::Local { cutoff } => {
            require_positive_finite(*cutoff.as_base(), "gamma::low_dose_cutoff")?
        }
    }

    let grid = *reference.grid();
    let [nx, ny, nz] = grid.dims();
    let spacing = grid.spacing();
    let zero = <T as NumericElement>::ZERO;

    let dta_base = *dta.as_base();
    let search_radius_base = *search_radius.as_base();
    let inv_dta_sq = (dta_base * dta_base).recip();
    let search_sq = search_radius_base * search_radius_base;

    // Per-axis neighborhood radius in voxels covering the search sphere.
    let radius_vox = |axis: usize| -> usize {
        (search_radius_base * spacing[axis].recip()).ceil().to_f64() as usize
    };
    let (rx, ry, rz) = (radius_vox(0), radius_vox(1), radius_vox(2));

    Ok(Volume::from_shape_fn(grid, |idx| {
        let (i, j, k) = (idx[0], idx[1], idx[2]);
        let world_r = grid.voxel_center(i, j, k);
        let dose_r = reference.get(i, j, k).expect("reference index within grid");

        // Per-point dose-difference denominator (global constant or local D_r);
        // low-dose points under local normalization are excluded (gamma 0).
        let delta_dose = match norm {
            Norm::Global(d) => dose_diff_criterion * *d.as_base(),
            Norm::Local { cutoff } => {
                if dose_r < *cutoff.as_base() {
                    return zero;
                }
                dose_diff_criterion * dose_r
            }
        };
        let inv_dd_sq = (delta_dose * delta_dose).recip();

        let mut best_sq = T::infinity();
        for ei in i.saturating_sub(rx)..=(i + rx).min(nx - 1) {
            for ej in j.saturating_sub(ry)..=(j + ry).min(ny - 1) {
                for ek in k.saturating_sub(rz)..=(k + rz).min(nz - 1) {
                    let world_e = grid.voxel_center(ei, ej, ek);
                    let dist_sq = (world_e - world_r).norm_squared();
                    if dist_sq > search_sq {
                        continue;
                    }
                    let dose_e = evaluated
                        .get(ei, ej, ek)
                        .expect("evaluated index within grid");
                    let dd = dose_e - dose_r;
                    let gamma_sq = dist_sq * inv_dta_sq + dd * dd * inv_dd_sq;
                    best_sq = best_sq.min_scalar(gamma_sq);
                }
            }
        }
        best_sq.sqrt()
    }))
}

/// Fraction of reference voxels above `dose_threshold` whose gamma is `≤ 1`.
///
/// `gamma` and `reference` must share a grid (both come from
/// [`gamma_index_3d`]'s inputs). If no voxel exceeds the threshold the pass rate
/// is defined as `1` (no failing points).
#[must_use]
pub fn gamma_pass_rate<T: Scalar>(
    gamma: &Volume<T>,
    reference: &Volume<T>,
    dose_threshold: AbsorbedDose<T>,
) -> T {
    let [nx, ny, nz] = reference.grid().dims();
    let one = <T as NumericElement>::ONE;
    let mut evaluated = 0usize;
    let mut passed = 0usize;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if reference.get(i, j, k).expect("index") < *dose_threshold.as_base() {
                    continue;
                }
                evaluated += 1;
                if gamma.get(i, j, k).expect("index") <= one {
                    passed += 1;
                }
            }
        }
    }
    if evaluated == 0 {
        return one;
    }
    T::from_f64(passed as f64) * T::from_f64(evaluated as f64).recip()
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    fn grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([5, 5, 5], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid")
    }

    fn distance(value: f64) -> Length<f64> {
        Length::from_base(value)
    }

    fn absorbed(value: f64) -> AbsorbedDose<f64> {
        AbsorbedDose::from_base(value)
    }

    #[test]
    fn identical_distributions_have_zero_gamma_and_full_pass() {
        let dose = Volume::from_shape_fn(grid(), |idx| 1.0 + idx[0] as f64);
        let gamma = gamma_index_3d(
            &dose,
            &dose,
            0.03,
            distance(2.0),
            absorbed(5.0),
            distance(4.0),
        )
        .expect("valid");
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    assert_relative_eq!(gamma.get(i, j, k).unwrap(), 0.0, epsilon = 1e-12);
                }
            }
        }
        assert_relative_eq!(
            gamma_pass_rate(&gamma, &dose, absorbed(0.0)),
            1.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn uniform_difference_scales_gamma_by_dose_ratio() {
        // Zero spatial gradient: the minimizing point is co-located, so the DTA
        // term is 0 and γ = |ΔD_actual| / ΔD_criterion. At exactly the criterion
        // γ = 1 (the analytical boundary value); at half the criterion γ = 0.5 and
        // the plan passes. (Pass/fail is asserted strictly off the γ=1 knife-edge,
        // which is rounding-sensitive by construction.)
        let norm = 10.0;
        let criterion_delta = 0.03 * norm;
        let reference = Volume::from_shape_fn(grid(), |_| norm);

        let at_criterion = Volume::from_shape_fn(grid(), |_| norm + criterion_delta);
        let gamma = gamma_index_3d(
            &reference,
            &at_criterion,
            0.03,
            distance(2.0),
            absorbed(norm),
            distance(4.0),
        )
        .expect("valid");
        assert_relative_eq!(gamma.get(2, 2, 2).unwrap(), 1.0, epsilon = 1e-12);

        let half = Volume::from_shape_fn(grid(), |_| norm + 0.5 * criterion_delta);
        let gamma_half = gamma_index_3d(
            &reference,
            &half,
            0.03,
            distance(2.0),
            absorbed(norm),
            distance(4.0),
        )
        .expect("valid");
        assert_relative_eq!(gamma_half.get(2, 2, 2).unwrap(), 0.5, epsilon = 1e-12);
        assert_relative_eq!(
            gamma_pass_rate(&gamma_half, &reference, absorbed(0.0)),
            1.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn difference_twice_the_criterion_fails_everywhere() {
        let norm = 10.0;
        let reference = Volume::from_shape_fn(grid(), |_| norm);
        let evaluated = Volume::from_shape_fn(grid(), |_| norm + 2.0 * 0.03 * norm);
        let gamma = gamma_index_3d(
            &reference,
            &evaluated,
            0.03,
            distance(2.0),
            absorbed(norm),
            distance(4.0),
        )
        .expect("valid");
        // γ = 2 everywhere (uniform, no closer match).
        assert_relative_eq!(gamma.get(2, 2, 2).unwrap(), 2.0, epsilon = 1e-12);
        assert_relative_eq!(
            gamma_pass_rate(&gamma, &reference, absorbed(0.0)),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn mismatched_grids_and_bad_criteria_error() {
        let a = Volume::from_shape_fn(grid(), |_| 1.0);
        let other = VoxelGrid::axis_aligned([4, 4, 4], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let b = Volume::from_shape_fn(other, |_| 1.0);
        assert!(gamma_index_3d(&a, &b, 0.03, distance(2.0), absorbed(5.0), distance(4.0)).is_err());
        assert!(gamma_index_3d(&a, &a, 0.0, distance(2.0), absorbed(5.0), distance(4.0)).is_err());
        assert!(
            gamma_index_3d(&a, &a, 0.03, distance(2.0), absorbed(-5.0), distance(4.0)).is_err()
        );
    }

    #[test]
    fn local_gamma_equals_global_for_uniform_dose() {
        // Uniform reference D=10, evaluated 3% high (0.3). Local ΔD = 0.03·10 = 0.3
        // = global ΔD with norm 10 → both give γ = 1 (co-located, no spatial term).
        let reference = Volume::from_shape_fn(grid(), |_| 10.0);
        let evaluated = Volume::from_shape_fn(grid(), |_| 10.3);
        let local = gamma_index_3d_local(
            &reference,
            &evaluated,
            0.03,
            distance(2.0),
            absorbed(1.0),
            distance(4.0),
        )
        .unwrap();
        assert_relative_eq!(local.get(2, 2, 2).unwrap(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn local_gamma_is_stricter_in_low_dose_than_global() {
        // Two spatial regions: i<2 → 10 Gy, i≥2 → 2 Gy; evaluated 3% high locally.
        // At an interior low-dose voxel, local γ = 1 (ΔD = 0.03·2) but global γ =
        // 0.2 (ΔD = 0.03·10) — local normalization is the stricter, correct metric.
        let reference = Volume::from_shape_fn(grid(), |idx| if idx[0] < 2 { 10.0 } else { 2.0 });
        let evaluated = Volume::from_shape_fn(grid(), |idx| {
            let d = if idx[0] < 2 { 10.0 } else { 2.0 };
            d * 1.03
        });
        let local = gamma_index_3d_local(
            &reference,
            &evaluated,
            0.03,
            distance(2.0),
            absorbed(0.5),
            distance(4.0),
        )
        .unwrap();
        let global = gamma_index_3d(
            &reference,
            &evaluated,
            0.03,
            distance(2.0),
            absorbed(10.0),
            distance(4.0),
        )
        .unwrap();
        // Interior low-dose voxel (i=3, neighbours i=2,4 also low).
        assert_relative_eq!(local.get(3, 2, 2).unwrap(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(global.get(3, 2, 2).unwrap(), 0.2, epsilon = 1e-12);
        assert!(local.get(3, 2, 2).unwrap() > global.get(3, 2, 2).unwrap());
    }

    #[test]
    fn local_gamma_excludes_points_below_the_low_dose_cutoff() {
        // Reference below the cutoff → excluded (gamma 0) even with a large error.
        let reference = Volume::from_shape_fn(grid(), |_| 0.5);
        let evaluated = Volume::from_shape_fn(grid(), |_| 5.0); // huge error
        let local = gamma_index_3d_local(
            &reference,
            &evaluated,
            0.03,
            distance(2.0),
            absorbed(1.0),
            distance(4.0),
        )
        .unwrap();
        assert_relative_eq!(local.get(2, 2, 2).unwrap(), 0.0, epsilon = 1e-15);
        // Bad cutoff is rejected.
        assert!(gamma_index_3d_local(
            &reference,
            &evaluated,
            0.03,
            distance(2.0),
            absorbed(-1.0),
            distance(4.0),
        )
        .is_err());
    }

    #[test]
    fn gamma_is_generic_over_scalar_f32() {
        let dose = Volume::from_shape_fn(
            VoxelGrid::<f32>::axis_aligned([3, 3, 3], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap(),
            |_| 1.0_f32,
        );
        let gamma = gamma_index_3d(
            &dose,
            &dose,
            0.03_f32,
            Length::from_base(2.0_f32),
            AbsorbedDose::from_base(5.0_f32),
            Length::from_base(4.0_f32),
        )
        .expect("valid");
        assert_relative_eq!(gamma.get(1, 1, 1).unwrap(), 0.0_f32, epsilon = 1e-6);
    }
}
