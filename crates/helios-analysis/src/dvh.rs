//! Cumulative dose-volume histogram (DVH).

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
    sorted: Vec<T>,
}

impl<T: Scalar> Dvh<T> {
    /// Build a DVH from every voxel of a dose volume.
    #[must_use]
    pub fn from_volume(dose: &Volume<T>) -> Self {
        Self::from_volume_masked(dose, |_| true)
    }

    /// Build a **structure-masked** DVH from the voxels of `dose` for which
    /// `include(idx)` is true — the per-structure (PTV / OAR) DVH clinical plan
    /// evaluation and DVH-agreement metrics operate on. [`from_volume`] is the
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
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if include([i, j, k]) {
                        sorted.push(dose.get(i, j, k).expect("index within grid"));
                    }
                }
            }
        }
        sorted.sort_by(|a, b| a.to_f64().total_cmp(&b.to_f64()));
        Self { sorted }
    }

    /// Number of voxels contributing to the histogram.
    #[must_use]
    pub fn count(&self) -> usize {
        self.sorted.len()
    }

    /// Minimum dose.
    #[must_use]
    pub fn min(&self) -> T {
        *self.sorted.first().expect("non-empty DVH")
    }

    /// Maximum dose.
    #[must_use]
    pub fn max(&self) -> T {
        *self.sorted.last().expect("non-empty DVH")
    }

    /// Mean dose over all sampled voxels.
    #[must_use]
    pub fn mean(&self) -> T {
        let sum = self
            .sorted
            .iter()
            .copied()
            .fold(<T as NumericElement>::ZERO, |acc, d| acc + d);
        sum * T::from_f64(self.sorted.len() as f64).recip()
    }

    /// Volume fraction (in `[0, 1]`) receiving **at least** `dose`.
    #[must_use]
    pub fn volume_fraction_at_dose(&self, dose: T) -> T {
        let at_least = self.sorted.iter().filter(|&&d| d >= dose).count();
        T::from_f64(at_least as f64) * T::from_f64(self.sorted.len() as f64).recip()
    }

    /// Near-rank dose `Dx`: the dose received by at least `fraction` of the
    /// volume (`fraction` in `[0, 1]`). `D_1.0` is the minimum dose, `D_0.0` the
    /// maximum.
    ///
    /// Nearest-rank (no interpolation): `k = ceil(fraction·n)` hottest voxels
    /// must meet the threshold, so `Dx` is the `k`-th largest dose.
    #[must_use]
    pub fn dose_at_volume_fraction(&self, fraction: T) -> T {
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
        if d50 <= <T as NumericElement>::ZERO {
            return <T as NumericElement>::ZERO;
        }
        (d2 - d98) * d50.recip()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
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
        assert_relative_eq!(dvh.min(), 2.5, epsilon = 1e-15);
        assert_relative_eq!(dvh.max(), 2.5, epsilon = 1e-15);
        assert_relative_eq!(dvh.mean(), 2.5, epsilon = 1e-15);
        // At/below 2.5 → full volume; above → none.
        assert_relative_eq!(dvh.volume_fraction_at_dose(2.5), 1.0, epsilon = 1e-15);
        assert_relative_eq!(dvh.volume_fraction_at_dose(2.5001), 0.0, epsilon = 1e-15);
        // Every Dx equals the uniform dose.
        assert_relative_eq!(dvh.dose_at_volume_fraction(0.5), 2.5, epsilon = 1e-15);
        assert_relative_eq!(dvh.dose_at_volume_fraction(1.0), 2.5, epsilon = 1e-15);
    }

    #[test]
    fn linear_ramp_has_known_mean_and_quantiles() {
        // 100 voxels with dose = index 0..99.
        let g = VoxelGrid::axis_aligned([100, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let dose = Volume::from_shape_fn(g, |idx| idx[0] as f64);
        let dvh = Dvh::from_volume(&dose);
        assert_relative_eq!(dvh.min(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(dvh.max(), 99.0, epsilon = 1e-12);
        assert_relative_eq!(dvh.mean(), 49.5, epsilon = 1e-12);
        // Half the volume (50 hottest voxels) receives ≥ dose 50 (values 50..99).
        assert_relative_eq!(dvh.dose_at_volume_fraction(0.5), 50.0, epsilon = 1e-12);
        // Exactly half of voxels have dose ≥ 50 (values 50..99 = 50 voxels).
        assert_relative_eq!(dvh.volume_fraction_at_dose(50.0), 0.5, epsilon = 1e-12);
    }

    #[test]
    fn masked_dvh_restricts_to_the_structure() {
        // 4×4×4 dose: a "target" (i<2) at 2.0 Gy, surrounding "OAR" (i≥2) at 8.0.
        let dose = Volume::from_shape_fn(grid([4, 4, 4]), |idx| if idx[0] < 2 { 2.0 } else { 8.0 });

        let target = Dvh::from_volume_masked(&dose, |idx| idx[0] < 2);
        assert_eq!(target.count(), 32); // half of 64 voxels
        assert_relative_eq!(target.mean(), 2.0, epsilon = 1e-15);
        assert_relative_eq!(target.max(), 2.0, epsilon = 1e-15);

        let oar = Dvh::from_volume_masked(&dose, |idx| idx[0] >= 2);
        assert_eq!(oar.count(), 32);
        assert_relative_eq!(oar.mean(), 8.0, epsilon = 1e-15);

        // Whole-volume mean (5.0) differs from either structure — masking matters.
        assert_relative_eq!(Dvh::from_volume(&dose).mean(), 5.0, epsilon = 1e-15);
    }

    #[test]
    fn single_voxel_mask_is_a_point_dvh() {
        let dose = Volume::from_shape_fn(grid([3, 3, 3]), |idx| (idx[0] + idx[1] + idx[2]) as f64);
        let point = Dvh::from_volume_masked(&dose, |idx| idx == [2, 2, 2]);
        assert_eq!(point.count(), 1);
        assert_relative_eq!(point.min(), 6.0, epsilon = 1e-15);
        assert_relative_eq!(point.dose_at_volume_fraction(1.0), 6.0, epsilon = 1e-15);
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
        assert_relative_eq!(dvh.mean(), 1.0_f32, epsilon = 1e-6);
    }
}
