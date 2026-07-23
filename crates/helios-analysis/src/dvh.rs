//! Cumulative dose-volume histogram (DVH).

use aequitas::systems::si::quantities::AbsorbedDose;
use asclepius::{
    response::radiation::{
        GeneralizedEquivalentUniformDose, LogisticControlProbability, LymanComplicationProbability,
    },
    BiologicalResponse, ResponseError, ResponseSlope, VolumeEffect,
};
use helios_domain::Volume;
use helios_math::{NumericElement, Scalar};

/// A cumulative dose-volume histogram built from a dose volume.
///
/// "Cumulative" means [`volume_fraction_at_dose`](Self::volume_fraction_at_dose)
/// reports the fraction of the sampled volume receiving **at least** a given
/// dose — the standard clinical DVH. Doses are retained (sorted ascending) so
/// quantile metrics (`Dx`, `Vx`) are exact under a nearest-rank convention.
#[derive(Debug, Clone)]
pub struct Dvh<T: Scalar> {
    /// Voxel doses, sorted ascending.
    sorted: Vec<AbsorbedDose<T>>,
    /// Whether the sample contains NaN values that are not ordered by `<`.
    contains_nan: bool,
}

impl<T: Scalar> Dvh<T> {
    /// Build a DVH from every voxel of a dose volume.
    #[must_use]
    pub fn from_volume(dose: &Volume<T>) -> Self {
        Self::from_volume_masked(dose, |_| true)
    }

    /// Build a **structure-masked** DVH from the voxels of `dose` for which
    /// `include(idx)` is true — the per-structure (PTV / OAR) DVH clinical plan
    /// evaluation and DVH-agreement metrics operate on. [`Self::from_volume`] is the
    /// whole-volume case (`include ≡ true`).
    ///
    /// The mask predicate is the segmentation contour (an ROI binary mask, e.g.
    /// from a `ritk` RT-struct rasterization) expressed as a voxel-index test.
    ///
    /// # Panics
    /// The quantile/statistic accessors require a non-empty histogram, so `include`
    /// must select at least one voxel.
    #[must_use]
    pub fn from_volume_masked<F>(dose: &Volume<T>, mut include: F) -> Self
    where
        F: FnMut([usize; 3]) -> bool,
    {
        let [nx, ny, nz] = dose.grid().dims();
        let mut sorted = Vec::new();
        let mut contains_nan = false;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if include([i, j, k]) {
                        let value = dose.get(i, j, k).expect("index within grid");
                        contains_nan |= value.to_f64().is_nan();
                        sorted.push(AbsorbedDose::from_base(value));
                    }
                }
            }
        }
        sorted.sort_by(|a, b| a.as_base().to_f64().total_cmp(&b.as_base().to_f64()));
        Self {
            sorted,
            contains_nan,
        }
    }

    /// Number of voxels contributing to the histogram.
    #[must_use]
    pub fn count(&self) -> usize {
        self.sorted.len()
    }

    /// Minimum dose as an Aequitas quantity.
    #[must_use]
    pub fn min(&self) -> AbsorbedDose<T> {
        *self.sorted.first().expect("non-empty DVH")
    }

    /// Maximum dose as an Aequitas quantity.
    #[must_use]
    pub fn max(&self) -> AbsorbedDose<T> {
        *self.sorted.last().expect("non-empty DVH")
    }

    /// Mean dose over all sampled voxels as an Aequitas quantity.
    #[must_use]
    pub fn mean(&self) -> AbsorbedDose<T> {
        let sum = self.sorted.iter().copied().fold(
            AbsorbedDose::from_base(<T as NumericElement>::ZERO),
            |acc, dose| acc + dose,
        );
        sum / T::from_f64(self.sorted.len() as f64)
    }

    /// Volume fraction (in `[0, 1]`) receiving **at least** `dose`.
    ///
    /// The threshold is an [`AbsorbedDose`] so a caller cannot pass an
    /// unrelated scalar physical quantity by accident.
    #[must_use]
    pub fn volume_fraction_at_dose(&self, dose: AbsorbedDose<T>) -> T {
        // NaN is unordered, so preserve the pre-indexed filter semantics for
        // invalid samples instead of treating a NaN suffix as qualifying. For
        // finite and infinite samples, the sorted invariant makes the lower
        // bound exact and reduces repeated queries from O(n) to O(log n) with
        // no allocation.
        let at_least = if dose.as_base().to_f64().is_nan() {
            0
        } else if self.contains_nan {
            self.sorted.iter().filter(|value| **value >= dose).count()
        } else {
            self.sorted.len() - self.sorted.partition_point(|value| *value < dose)
        };
        T::from_f64(at_least as f64) * T::from_f64(self.sorted.len() as f64).recip()
    }

    /// Near-rank dose `Dx`: the dose received by at least `fraction` of the
    /// volume (`fraction` in `[0, 1]`). `D_1.0` is the minimum dose, `D_0.0` the
    /// maximum. The returned dose retains its Aequitas dimension.
    ///
    /// Nearest-rank (no interpolation): `k = ceil(fraction·n)` hottest voxels
    /// must meet the threshold, so `Dx` is the `k`-th largest dose.
    #[must_use]
    pub fn dose_at_volume_fraction(&self, fraction: T) -> AbsorbedDose<T> {
        let n = self.sorted.len();
        let frac = fraction.to_f64().clamp(0.0, 1.0);
        // k hottest voxels; k in [1, n]. k=0 (fraction 0) → hottest voxel.
        let k = (frac * n as f64).ceil() as usize;
        let k = k.clamp(1, n);
        self.sorted[n - k]
    }

    /// ICRU-83 dose **homogeneity index** `HI = (D₂ − D₉₈) / D₅₀` over the sampled
    /// (usually target) volume. Lower is more homogeneous; a perfectly uniform
    /// dose gives `0`. Returns `0` when `D₅₀` is zero (no dose to normalize by).
    #[must_use]
    pub fn homogeneity_index(&self) -> T {
        let d2 = self.dose_at_volume_fraction(T::from_f64(0.02));
        let d98 = self.dose_at_volume_fraction(T::from_f64(0.98));
        let d50 = self.dose_at_volume_fraction(T::from_f64(0.5));
        if *d50.as_base() <= <T as NumericElement>::ZERO {
            return <T as NumericElement>::ZERO;
        }
        ((d2 - d98) / d50).into_base()
    }

    /// The structure's dose sample (ascending-sorted), borrowed zero-copy.
    ///
    /// The same voxel doses the histogram summarizes; the radiobiology metrics
    /// below operate on this sample without re-scanning the dose volume.
    #[must_use]
    pub fn dose_sample(&self) -> &[AbsorbedDose<T>] {
        &self.sorted
    }

    /// Generalized equivalent uniform dose (gEUD) of this structure's dose, with
    /// volume-effect parameter `a`. The canonical Asclepius law borrows the
    /// already-sampled Aequitas dose quantities, with no allocation or volume
    /// re-scan.
    ///
    /// # Errors
    ///
    /// Returns [`ResponseError`] when `a` is zero or non-finite, the sample is
    /// empty, a dose is negative or non-finite, or a negative exponent observes
    /// zero dose.
    pub fn generalized_eud(&self, a: T) -> Result<AbsorbedDose<T>, ResponseError<T>> {
        let volume_effect = VolumeEffect::new(a).map_err(ResponseError::from)?;
        GeneralizedEquivalentUniformDose::new(volume_effect).evaluate(&self.sorted)
    }

    /// Niemierko logistic tumour control probability of this structure's dose:
    /// the Asclepius logistic law evaluated at this structure's gEUD.
    ///
    /// # Errors
    ///
    /// Returns [`ResponseError`] when any model parameter or dose observation
    /// violates the law's mathematical domain.
    pub fn tcp_logistic(
        &self,
        a: T,
        tcd50: AbsorbedDose<T>,
        gamma50: T,
    ) -> Result<T, ResponseError<T>> {
        let dose = self.generalized_eud(a)?;
        let slope = ResponseSlope::new(gamma50).map_err(ResponseError::from)?;
        let model = LogisticControlProbability::new(tcd50, slope).map_err(ResponseError::from)?;
        model.evaluate(dose).map(asclepius::Probability::get)
    }

    /// Lyman–Kutcher–Burman normal-tissue complication probability of this
    /// structure's dose, evaluated through the canonical Asclepius law.
    ///
    /// # Errors
    ///
    /// Returns [`ResponseError`] when any model parameter or dose observation
    /// violates the law's mathematical domain.
    pub fn ntcp_lkb(&self, a: T, td50: AbsorbedDose<T>, m: T) -> Result<T, ResponseError<T>> {
        let dose = self.generalized_eud(a)?;
        let slope = ResponseSlope::new(m).map_err(ResponseError::from)?;
        let model = LymanComplicationProbability::new(td50, slope).map_err(ResponseError::from)?;
        model.evaluate(dose).map(asclepius::Probability::get)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    fn grid(dims: [usize; 3]) -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned(dims, [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0)).expect("grid")
    }

    #[test]
    fn uniform_dose_is_a_step_histogram() {
        let dose = Volume::from_shape_fn(grid([4, 4, 4]), |_| 2.5);
        let dvh = Dvh::from_volume(&dose);
        assert_eq!(dvh.count(), 64);
        assert_relative_eq!(dvh.min().into_base(), 2.5, epsilon = 1e-15);
        assert_relative_eq!(dvh.max().into_base(), 2.5, epsilon = 1e-15);
        assert_relative_eq!(dvh.mean().into_base(), 2.5, epsilon = 1e-15);
        // At/below 2.5 → full volume; above → none.
        assert_relative_eq!(
            dvh.volume_fraction_at_dose(AbsorbedDose::from_base(2.5)),
            1.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            dvh.volume_fraction_at_dose(AbsorbedDose::from_base(2.5001)),
            0.0,
            epsilon = 1e-15
        );
        // Every Dx equals the uniform dose.
        assert_relative_eq!(
            dvh.dose_at_volume_fraction(0.5).into_base(),
            2.5,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            dvh.dose_at_volume_fraction(1.0).into_base(),
            2.5,
            epsilon = 1e-15
        );
    }

    #[test]
    fn linear_ramp_has_known_mean_and_quantiles() {
        // 100 voxels with dose = index 0..99.
        let g = VoxelGrid::axis_aligned([100, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let dose = Volume::from_shape_fn(g, |idx| idx[0] as f64);
        let dvh = Dvh::from_volume(&dose);
        assert_relative_eq!(dvh.min().into_base(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(dvh.max().into_base(), 99.0, epsilon = 1e-12);
        assert_relative_eq!(dvh.mean().into_base(), 49.5, epsilon = 1e-12);
        // Half the volume (50 hottest voxels) receives ≥ dose 50 (values 50..99).
        assert_relative_eq!(
            dvh.dose_at_volume_fraction(0.5).into_base(),
            50.0,
            epsilon = 1e-12
        );
        // Exactly half of voxels have dose ≥ 50 (values 50..99 = 50 voxels).
        assert_relative_eq!(
            dvh.volume_fraction_at_dose(AbsorbedDose::from_base(50.0)),
            0.5,
            epsilon = 1e-12
        );
    }

    #[test]
    fn masked_dvh_restricts_to_the_structure() {
        // 4×4×4 dose: a "target" (i<2) at 2.0 Gy, surrounding "OAR" (i≥2) at 8.0.
        let dose = Volume::from_shape_fn(grid([4, 4, 4]), |idx| if idx[0] < 2 { 2.0 } else { 8.0 });

        let target = Dvh::from_volume_masked(&dose, |idx| idx[0] < 2);
        assert_eq!(target.count(), 32); // half of 64 voxels
        assert_relative_eq!(target.mean().into_base(), 2.0, epsilon = 1e-15);
        assert_relative_eq!(target.max().into_base(), 2.0, epsilon = 1e-15);

        let oar = Dvh::from_volume_masked(&dose, |idx| idx[0] >= 2);
        assert_eq!(oar.count(), 32);
        assert_relative_eq!(oar.mean().into_base(), 8.0, epsilon = 1e-15);

        // Whole-volume mean (5.0) differs from either structure — masking matters.
        assert_relative_eq!(
            Dvh::from_volume(&dose).mean().into_base(),
            5.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn single_voxel_mask_is_a_point_dvh() {
        let dose = Volume::from_shape_fn(grid([3, 3, 3]), |idx| (idx[0] + idx[1] + idx[2]) as f64);
        let point = Dvh::from_volume_masked(&dose, |idx| idx == [2, 2, 2]);
        assert_eq!(point.count(), 1);
        assert_relative_eq!(point.min().into_base(), 6.0, epsilon = 1e-15);
        assert_relative_eq!(
            point.dose_at_volume_fraction(1.0).into_base(),
            6.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn homogeneity_index_is_zero_for_uniform_and_known_for_a_ramp() {
        // Uniform dose → every Dx equal → HI = 0.
        let uniform = Dvh::from_volume(&Volume::from_shape_fn(grid([4, 4, 4]), |_| 3.0));
        assert_relative_eq!(uniform.homogeneity_index(), 0.0, epsilon = 1e-15);

        // Ramp 0..99 over 100 voxels: D2 = 98, D98 = 2, D50 = 50 → HI = 96/50 = 1.92.
        let g = VoxelGrid::axis_aligned([100, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let ramp = Dvh::from_volume(&Volume::from_shape_fn(g, |idx| idx[0] as f64));
        assert_relative_eq!(ramp.homogeneity_index(), 96.0 / 50.0, epsilon = 1e-12);
    }

    #[test]
    fn dvh_is_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([2, 2, 2], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .expect("grid");
        let dose = Volume::from_shape_fn(g, |_| 1.0_f32);
        let dvh = Dvh::from_volume(&dose);
        assert_relative_eq!(dvh.mean().into_base(), 1.0_f32, epsilon = 1e-6);
    }

    #[test]
    fn dvh_geud_reuses_the_sample_and_matches_the_canonical_law() {
        // Heterogeneous dose; the DVH method must equal direct Asclepius
        // evaluation over the identical borrowed sample.
        let dose = Volume::from_shape_fn(grid([4, 4, 4]), |idx| {
            1.0 + idx[0] as f64 + 0.5 * idx[1] as f64
        });
        let dvh = Dvh::from_volume(&dose);
        assert_eq!(dvh.dose_sample().len(), 64);
        for a in [1.0, 2.0, -3.0, 8.0] {
            let model = GeneralizedEquivalentUniformDose::new(
                VolumeEffect::new(a).expect("finite non-zero exponent"),
            );
            let direct = model
                .evaluate(dvh.dose_sample())
                .expect("valid dose observation");
            assert_relative_eq!(
                dvh.generalized_eud(a).expect("valid response").into_base(),
                direct.into_base(),
                max_relative = 1e-13
            );
        }
    }

    #[test]
    fn threshold_query_preserves_boundary_and_nan_semantics() {
        let g = VoxelGrid::axis_aligned([4, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        let finite = Volume::from_shape_fn(g, |[i, _, _]| i as f64);
        let finite_dvh = Dvh::from_volume(&finite);
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(-1.0)),
            1.0
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(0.0)),
            1.0
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(3.0)),
            0.25
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(4.0)),
            0.0
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(f64::NAN)),
            0.0
        );

        let with_nan = Volume::from_shape_fn(
            VoxelGrid::axis_aligned([3, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .expect("grid"),
            |[i, _, _]| [1.0, f64::NAN, 3.0][i],
        );
        let nan_dvh = Dvh::from_volume(&with_nan);
        assert_eq!(
            nan_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(1.0)),
            2.0 / 3.0
        );
        assert_eq!(
            nan_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(2.0)),
            1.0 / 3.0
        );
        assert!(matches!(
            nan_dvh.generalized_eud(1.0),
            Err(ResponseError::InvalidObservation { index: 2, .. })
        ));
    }

    #[test]
    fn outcome_parameters_preserve_asclepius_validation_errors() {
        let dvh = Dvh::from_volume(&Volume::from_shape_fn(grid([2, 1, 1]), |_| 1.0));
        assert!(matches!(
            dvh.generalized_eud(0.0),
            Err(ResponseError::InvalidValue(source))
                if source.kind() == asclepius::ValueKind::VolumeEffect
                    && source.constraint() == asclepius::ValueConstraint::FiniteNonZero
        ));
        assert!(matches!(
            dvh.tcp_logistic(1.0, AbsorbedDose::from_base(1.0), -0.5),
            Err(ResponseError::InvalidValue(source))
                if source.kind() == asclepius::ValueKind::ResponseSlope
                    && source.constraint() == asclepius::ValueConstraint::FinitePositive
        ));
    }

    #[test]
    fn uniform_structure_outcome_equals_the_pointwise_model() {
        // A uniform-dose structure has gEUD = that dose for any a, so its TCP/NTCP
        // reduce to the pointwise outcome models at that dose.
        let d = 62.0;
        let dvh = Dvh::from_volume(&Volume::from_shape_fn(grid([3, 3, 3]), |_| d));
        assert_relative_eq!(
            dvh.generalized_eud(-10.0)
                .expect("valid response")
                .into_base(),
            d,
            max_relative = 1e-12
        );
        // NTCP with TD50 = d ⇒ t = 0 ⇒ 0.5; TCP with TCD50 = d ⇒ 0.5.
        assert_relative_eq!(
            dvh.ntcp_lkb(1.0, AbsorbedDose::from_base(d), 0.2)
                .expect("valid response"),
            0.5,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            dvh.tcp_logistic(1.0, AbsorbedDose::from_base(d), 2.0)
                .expect("valid response"),
            0.5,
            epsilon = 1e-12
        );
        let geud = dvh.generalized_eud(2.0).expect("valid response");
        let ntcp_model = LymanComplicationProbability::new(
            AbsorbedDose::from_base(50.0),
            ResponseSlope::new(0.2).expect("positive slope"),
        )
        .expect("positive midpoint");
        let tcp_model = LogisticControlProbability::new(
            AbsorbedDose::from_base(55.0),
            ResponseSlope::new(2.0).expect("positive slope"),
        )
        .expect("positive midpoint");
        assert_relative_eq!(
            dvh.ntcp_lkb(2.0, AbsorbedDose::from_base(50.0), 0.2)
                .expect("valid response"),
            ntcp_model.evaluate(geud).expect("valid response").get(),
            epsilon = 1e-14
        );
        assert_relative_eq!(
            dvh.tcp_logistic(2.0, AbsorbedDose::from_base(55.0), 2.0)
                .expect("valid response"),
            tcp_model.evaluate(geud).expect("valid response").get(),
            epsilon = 1e-14
        );
    }

    #[test]
    fn masked_structure_outcome_reflects_only_masked_voxels() {
        // Two half-slabs at 20 and 80 Gy; masking the hot half gives gEUD ≈ 80,
        // higher NTCP than masking the cold half — the per-structure evaluation.
        let dose =
            Volume::from_shape_fn(grid([4, 4, 4]), |idx| if idx[0] < 2 { 20.0 } else { 80.0 });
        let hot = Dvh::from_volume_masked(&dose, |idx| idx[0] >= 2);
        let cold = Dvh::from_volume_masked(&dose, |idx| idx[0] < 2);
        assert_relative_eq!(
            hot.generalized_eud(1.0)
                .expect("valid response")
                .into_base(),
            80.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            cold.generalized_eud(1.0)
                .expect("valid response")
                .into_base(),
            20.0,
            epsilon = 1e-12
        );
        assert!(
            hot.ntcp_lkb(1.0, AbsorbedDose::from_base(50.0), 0.2)
                .expect("valid response")
                > cold
                    .ntcp_lkb(1.0, AbsorbedDose::from_base(50.0), 0.2)
                    .expect("valid response")
        );
    }
}
