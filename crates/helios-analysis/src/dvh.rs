//! Cumulative dose-volume histogram (DVH).

use aequitas::systems::si::quantities::{AbsorbedDose, Dimensionless};
use asclepius::{
    response::radiation::{
        GeneralizedEquivalentUniformDose, LogisticControlProbability, LymanComplicationProbability,
    },
    BiologicalResponse, Gamma50, LymanSlope, ResponseError, VolumeEffect,
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
    /// dimensionless volume-effect parameter `a`. The canonical Asclepius law borrows the
    /// already-sampled Aequitas dose quantities, with no allocation or volume
    /// re-scan.
    ///
    /// # Errors
    ///
    /// Returns [`ResponseError`] when `a` is zero or non-finite, the sample is
    /// empty, a dose is negative or non-finite, or a negative exponent observes
    /// zero dose.
    pub fn generalized_eud(
        &self,
        a: Dimensionless<T>,
    ) -> Result<AbsorbedDose<T>, ResponseError<T>> {
        let volume_effect = VolumeEffect::new(*a.as_base()).map_err(ResponseError::from)?;
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
        a: Dimensionless<T>,
        tcd50: AbsorbedDose<T>,
        gamma50: T,
    ) -> Result<T, ResponseError<T>> {
        let dose = self.generalized_eud(a)?;
        let gamma50 = Gamma50::new(gamma50).map_err(ResponseError::from)?;
        let model = LogisticControlProbability::new(tcd50, gamma50).map_err(ResponseError::from)?;
        model.evaluate(dose).map(asclepius::Probability::get)
    }

    /// Lyman–Kutcher–Burman normal-tissue complication probability of this
    /// structure's dose, evaluated through the canonical Asclepius law.
    ///
    /// # Errors
    ///
    /// Returns [`ResponseError`] when any model parameter or dose observation
    /// violates the law's mathematical domain.
    pub fn ntcp_lkb(
        &self,
        a: Dimensionless<T>,
        td50: AbsorbedDose<T>,
        m: T,
    ) -> Result<T, ResponseError<T>> {
        let dose = self.generalized_eud(a)?;
        let m = LymanSlope::new(m).map_err(ResponseError::from)?;
        let model = LymanComplicationProbability::new(td50, m).map_err(ResponseError::from)?;
        model.evaluate(dose).map(asclepius::Probability::get)
    }
}

#[cfg(test)]
mod tests;
