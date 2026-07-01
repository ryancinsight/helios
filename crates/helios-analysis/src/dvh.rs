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
        let [nx, ny, nz] = dose.grid().dims();
        let mut sorted = Vec::with_capacity(nx * ny * nz);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    sorted.push(dose.get(i, j, k).expect("index within grid"));
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
    fn dvh_is_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([2, 2, 2], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .expect("grid");
        let dose = Volume::from_shape_fn(g, |_| 1.0_f32);
        let dvh = Dvh::from_volume(&dose);
        assert_relative_eq!(dvh.mean(), 1.0_f32, epsilon = 1e-6);
    }
}
