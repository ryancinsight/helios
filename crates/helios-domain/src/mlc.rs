//! Binary multi-leaf collimator (MLC) modulation with leakage and
//! tongue-and-groove.
//!
//! TomoTherapy modulates the fan beam with a **binary** MLC: each leaf is open or
//! closed, and modulation is achieved by the *fraction* of each projection a leaf
//! is open (leaf-open-time, LOT), stored as a [`LeafOpenTimeSinogram`]. The
//! delivered fluence departs from the nominal open fraction through two effects
//! modelled here ([`MlcModel`]):
//!
//! - **Leakage / transmission**: a closed leaf still transmits a small fraction
//!   `τ` (inter-/intra-leaf leakage, ~0.5–1 % for TomoTherapy), so the effective
//!   fluence is `open + (1−open)·τ`.
//! - **Tongue-and-groove**: the overlapping tongue/groove at a leaf's sides
//!   underdoses the edge region where a neighbour is more closed, reducing the
//!   effective fluence by a first-order edge term.

use helios_core::HeliosError;
use helios_math::{NumericElement, Scalar};

/// Per-projection, per-leaf open-time fractions in `[0, 1]`, row-major
/// `[projection][leaf]`.
#[derive(Debug, Clone, PartialEq)]
pub struct LeafOpenTimeSinogram<T: Scalar> {
    projections: usize,
    leaves: usize,
    fractions: Vec<T>,
}

impl<T: Scalar> LeafOpenTimeSinogram<T> {
    /// Build from a flat `[projection][leaf]` vector of open-time fractions.
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if the length does not equal
    /// `projections·leaves`, or any fraction is non-finite or outside `[0, 1]`.
    pub fn from_fractions(
        projections: usize,
        leaves: usize,
        fractions: Vec<T>,
    ) -> Result<Self, HeliosError> {
        if fractions.len() != projections * leaves {
            return Err(HeliosError::InvalidDomainValue {
                field: "LeafOpenTimeSinogram::len",
                value: fractions.len() as f64,
                reason: "fraction count must equal projections·leaves",
            });
        }
        let (zero, one) = (<T as NumericElement>::ZERO, <T as NumericElement>::ONE);
        for &f in &fractions {
            if !f.is_finite() || f < zero || f > one {
                return Err(HeliosError::InvalidDomainValue {
                    field: "LeafOpenTimeSinogram::fraction",
                    value: f.to_f64(),
                    reason: "leaf open-time fraction must be finite in [0, 1]",
                });
            }
        }
        Ok(Self {
            projections,
            leaves,
            fractions,
        })
    }

    /// `(projections, leaves)`.
    #[must_use]
    pub fn dims(&self) -> (usize, usize) {
        (self.projections, self.leaves)
    }

    /// Open-time fraction of `leaf` at `projection`.
    #[must_use]
    pub fn get(&self, projection: usize, leaf: usize) -> T {
        self.fractions[projection * self.leaves + leaf]
    }
}

/// Binary-MLC fluence model: leakage transmission and tongue-and-groove.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MlcModel<T: Scalar> {
    leakage: T,
    tongue_and_groove: T,
}

impl<T: Scalar> MlcModel<T> {
    /// Construct from a leakage/transmission fraction (`[0, 1)`) and a
    /// tongue-and-groove edge fraction (`[0, 1]`).
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if either is non-finite or out
    /// of range.
    pub fn new(leakage: T, tongue_and_groove: T) -> Result<Self, HeliosError> {
        let (zero, one) = (<T as NumericElement>::ZERO, <T as NumericElement>::ONE);
        if !leakage.is_finite() || leakage < zero || leakage >= one {
            return Err(HeliosError::InvalidDomainValue {
                field: "MlcModel::leakage",
                value: leakage.to_f64(),
                reason: "leakage must be finite in [0, 1)",
            });
        }
        if !tongue_and_groove.is_finite() || tongue_and_groove < zero || tongue_and_groove > one {
            return Err(HeliosError::InvalidDomainValue {
                field: "MlcModel::tongue_and_groove",
                value: tongue_and_groove.to_f64(),
                reason: "tongue-and-groove fraction must be finite in [0, 1]",
            });
        }
        Ok(Self {
            leakage,
            tongue_and_groove,
        })
    }

    /// Leakage/transmission fraction.
    #[must_use]
    pub fn leakage(&self) -> T {
        self.leakage
    }

    /// Tongue-and-groove edge fraction.
    #[must_use]
    pub fn tongue_and_groove(&self) -> T {
        self.tongue_and_groove
    }

    /// Effective transmitted fluence for a leaf with open-time `open` whose
    /// neighbours are open `left`/`right` (in `[0, 1]`).
    ///
    /// `leakage-adjusted transmission − tongue-and-groove edge loss`, clamped to
    /// `[0, 1]`. The edge loss is `τ_tg · ½·(max(0, open−left) + max(0, open−right))`
    /// — the leaf loses fluence at each side where its neighbour is more closed.
    #[must_use]
    pub fn effective_fluence(&self, open: T, left: T, right: T) -> T {
        let (zero, one) = (<T as NumericElement>::ZERO, <T as NumericElement>::ONE);
        let base = open + (one - open) * self.leakage;
        let half = T::from_f64(0.5);
        let left_step = (open - left).max_scalar(zero);
        let right_step = (open - right).max_scalar(zero);
        let tg_loss = self.tongue_and_groove * half * (left_step + right_step);
        (base - tg_loss).max_scalar(zero).min_scalar(one)
    }

    /// Apply the model to a leaf-open-time sinogram, returning the effective
    /// fluence per `[projection][leaf]` (row-major). Outermost leaves have a
    /// single neighbour (the other side sees no tongue-and-groove partner).
    #[must_use]
    pub fn effective_fluence_sinogram(&self, lot: &LeafOpenTimeSinogram<T>) -> Vec<T> {
        let (projections, leaves) = lot.dims();
        let mut out = Vec::with_capacity(projections * leaves);
        for p in 0..projections {
            for l in 0..leaves {
                let open = lot.get(p, l);
                // Missing neighbour → same open (no T&G partner on that side).
                let left = if l > 0 { lot.get(p, l - 1) } else { open };
                let right = if l + 1 < leaves {
                    lot.get(p, l + 1)
                } else {
                    open
                };
                out.push(self.effective_fluence(open, left, right));
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn sinogram_validates_length_and_range() {
        assert!(LeafOpenTimeSinogram::from_fractions(2, 3, vec![0.0; 6]).is_ok());
        assert!(LeafOpenTimeSinogram::from_fractions(2, 3, vec![0.0; 5]).is_err());
        assert!(LeafOpenTimeSinogram::from_fractions(1, 1, vec![1.5]).is_err());
        assert!(LeafOpenTimeSinogram::from_fractions(1, 1, vec![-0.1]).is_err());
    }

    #[test]
    fn model_rejects_out_of_range() {
        assert!(MlcModel::new(1.0_f64, 0.1).is_err()); // leakage must be < 1
        assert!(MlcModel::new(-0.1_f64, 0.1).is_err());
        assert!(MlcModel::new(0.01_f64, 1.5).is_err());
        assert!(MlcModel::new(0.01_f64, 0.1).is_ok());
    }

    #[test]
    fn leakage_sets_closed_and_open_extremes() {
        let m = MlcModel::new(0.01_f64, 0.0).unwrap();
        // Fully closed leaf still transmits the leakage fraction.
        assert_relative_eq!(m.effective_fluence(0.0, 0.0, 0.0), 0.01, epsilon = 1e-15);
        // Fully open leaf transmits fully.
        assert_relative_eq!(m.effective_fluence(1.0, 1.0, 1.0), 1.0, epsilon = 1e-15);
        // Half-open: 0.5 + 0.5·0.01 = 0.505.
        assert_relative_eq!(m.effective_fluence(0.5, 0.5, 0.5), 0.505, epsilon = 1e-15);
    }

    #[test]
    fn tongue_and_groove_underdoses_isolated_open_leaves() {
        let m = MlcModel::new(0.0_f64, 0.1).unwrap();
        // Open leaf between open neighbours: no T&G loss.
        assert_relative_eq!(m.effective_fluence(1.0, 1.0, 1.0), 1.0, epsilon = 1e-15);
        // Open leaf between two closed neighbours: loss = 0.1·½·(1+1) = 0.1.
        assert_relative_eq!(m.effective_fluence(1.0, 0.0, 0.0), 0.9, epsilon = 1e-15);
        // One closed side only: loss = 0.1·½·(0+1) = 0.05.
        assert_relative_eq!(m.effective_fluence(1.0, 1.0, 0.0), 0.95, epsilon = 1e-15);
        // A more-open neighbour never adds loss (max(0, ·)).
        assert_relative_eq!(m.effective_fluence(0.3, 1.0, 1.0), 0.3, epsilon = 1e-15);
    }

    #[test]
    fn effective_fluence_stays_in_unit_interval() {
        let m = MlcModel::new(0.02_f64, 0.5).unwrap();
        for &(o, l, r) in &[(1.0, 0.0, 0.0), (0.0, 1.0, 1.0), (0.7, 0.1, 0.9)] {
            let f = m.effective_fluence(o, l, r);
            assert!((0.0..=1.0).contains(&f), "f={f}");
        }
    }

    #[test]
    fn sinogram_application_uses_neighbours() {
        // Row: [1, 0, 1] with T&G 0.1, no leakage.
        // leaf0: open=1, left=self(1), right=0 → loss 0.1·½·(0+1)=0.05 → 0.95
        // leaf1: open=0 → base 0, no loss → 0
        // leaf2: open=1, left=0, right=self(1) → loss 0.05 → 0.95
        let lot = LeafOpenTimeSinogram::from_fractions(1, 3, vec![1.0, 0.0, 1.0]).unwrap();
        let m = MlcModel::new(0.0_f64, 0.1).unwrap();
        let eff = m.effective_fluence_sinogram(&lot);
        assert_relative_eq!(eff[0], 0.95, epsilon = 1e-15);
        assert_relative_eq!(eff[1], 0.0, epsilon = 1e-15);
        assert_relative_eq!(eff[2], 0.95, epsilon = 1e-15);
    }

    #[test]
    fn mlc_is_generic_over_scalar_f32() {
        let m = MlcModel::new(0.01_f32, 0.1).unwrap();
        assert_relative_eq!(m.effective_fluence(0.0_f32, 0.0, 0.0), 0.01, epsilon = 1e-6);
        assert_relative_eq!(m.effective_fluence(1.0_f32, 0.0, 0.0), 0.9, epsilon = 1e-6);
    }
}
