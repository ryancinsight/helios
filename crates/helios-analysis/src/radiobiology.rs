//! Radiobiological plan-evaluation metrics: generalized equivalent uniform dose
//! (gEUD) and the TCP / NTCP outcome models built on it.
//!
//! These summarize a dose distribution into a single dose-response prediction,
//! the complement to the spatial DVH / gamma metrics: [`generalized_eud`]
//! collapses a dose sample to one biologically-weighted dose, and
//! [`tcp_logistic`] / [`ntcp_lkb`] map that gEUD to a tumour-control /
//! normal-tissue-complication probability. All are generic over the [`Scalar`]
//! seam and are pure functions (no allocation), so they compose over any ROI's
//! dose sample.

use helios_math::{NumericElement, Scalar};

/// Niemierko generalized equivalent uniform dose (gEUD) of a dose sample:
///
/// ```text
/// gEUD = ( (1/N) Σ_i D_i^a )^(1/a)
/// ```
///
/// The volume-effect parameter `a` tunes the tissue response: `a = 1` gives the
/// mean dose, `a → +∞` the maximum (serial organs / hot-spot control), and
/// `a → −∞` the minimum (parallel targets / cold-spot coverage). gEUD is the
/// Niemierko power mean, so it is bounded in `[min, max]` and monotonically
/// increasing in `a`; a uniform dose yields that dose for any `a`.
///
/// Doses must be non-negative (a non-integer `a` needs `D_i > 0` for a finite
/// result), and `a ≠ 0`. An empty sample returns zero.
///
/// # Panics
/// Debug-asserts `a ≠ 0`; in release `a = 0` yields a non-finite result.
#[must_use]
pub fn generalized_eud<T: Scalar>(doses: &[T], a: T) -> T {
    let zero = <T as NumericElement>::ZERO;
    debug_assert!(a != zero, "gEUD volume parameter a must be non-zero");
    if doses.is_empty() {
        return zero;
    }
    let inv_n = T::from_f64(doses.len() as f64).recip();
    let mean_pow = doses.iter().fold(zero, |acc, &d| acc + d.powf(a)) * inv_n;
    mean_pow.powf(a.recip())
}

/// Niemierko logistic tumour control probability (TCP) as a function of gEUD:
///
/// ```text
/// TCP = 1 / ( 1 + (TCD50 / gEUD)^(4·γ50) )
/// ```
///
/// `tcd50` is the gEUD giving 50 % control and `gamma50` the normalized
/// dose-response slope at that point. TCP is bounded in `[0, 1]`, equals `0.5`
/// at `gEUD = tcd50` for any `gamma50`, and increases monotonically with gEUD
/// (→ 1 as gEUD → ∞, → 0 as gEUD → 0). `gEUD` and `tcd50` must be positive.
#[must_use]
pub fn tcp_logistic<T: Scalar>(geud: T, tcd50: T, gamma50: T) -> T {
    let one = <T as NumericElement>::ONE;
    let four = T::from_f64(4.0);
    let ratio = (tcd50 * geud.recip()).powf(four * gamma50);
    (one + ratio).recip()
}

/// Lyman–Kutcher–Burman (LKB) normal-tissue complication probability (NTCP) as a
/// function of gEUD:
///
/// ```text
/// t    = (gEUD − TD50) / (m · TD50)
/// NTCP = Φ(t) = ½ · erfc(−t / √2)
/// ```
///
/// where `Φ` is the standard normal CDF. `td50` is the gEUD giving 50 %
/// complication and `m` the slope parameter (smaller `m` ⇒ steeper). NTCP is
/// bounded in `[0, 1]`, equals `0.5` at `gEUD = td50`, and increases
/// monotonically with gEUD. `td50` and `m` must be positive.
#[must_use]
pub fn ntcp_lkb<T: Scalar>(geud: T, td50: T, m: T) -> T {
    let half = T::from_f64(0.5);
    let inv_sqrt2 = T::from_f64(core::f64::consts::FRAC_1_SQRT_2);
    let t = (geud - td50) * (m * td50).recip();
    half * (-t * inv_sqrt2).erfc()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn generalized_eud_recovers_known_limits() {
        // a = 1 is the arithmetic mean.
        assert_relative_eq!(generalized_eud(&[1.0, 2.0, 3.0], 1.0), 2.0, epsilon = 1e-13);
        // A uniform dose is its own gEUD for any a.
        for a in [2.0, -4.0, 8.0, 0.5] {
            assert_relative_eq!(
                generalized_eud(&[5.0, 5.0, 5.0], a),
                5.0,
                max_relative = 1e-12
            );
        }
        // gEUD is the power mean: bounded in [min, max], increasing in a. It nears
        // the extremes only slowly ((1/N)^{1/a}) and large |a| overflows, so assert
        // the rigorous bounds + ordering, not tight equality with max/min.
        let d = [1.0, 2.0, 10.0]; // min 1, mean ≈ 4.333, max 10
        let (lo, hi) = (generalized_eud(&d, -40.0), generalized_eud(&d, 40.0));
        assert!(
            (1.0..=10.0).contains(&lo),
            "gEUD(a=−40) {lo} out of [min,max]"
        );
        assert!(
            (1.0..=10.0).contains(&hi),
            "gEUD(a=+40) {hi} out of [min,max]"
        );
        let mean = generalized_eud(&d, 1.0);
        assert!(
            lo < mean && mean < hi,
            "gEUD must increase with a: {lo} < {mean} < {hi}"
        );
        assert!(lo < 1.5 && hi > 9.0, "gEUD should near min/max: {lo}, {hi}");
    }

    #[test]
    fn tcp_is_bounded_half_at_tcd50_and_monotone() {
        let (tcd50, gamma50) = (60.0, 2.0);
        assert_relative_eq!(tcp_logistic(tcd50, tcd50, gamma50), 0.5, epsilon = 1e-12);
        let (low, high) = (
            tcp_logistic(40.0, tcd50, gamma50),
            tcp_logistic(80.0, tcd50, gamma50),
        );
        assert!(
            0.0 < low && low < 0.5,
            "TCP below TCD50 must be in (0,0.5): {low}"
        );
        assert!(
            0.5 < high && high < 1.0,
            "TCP above TCD50 must be in (0.5,1): {high}"
        );
        assert!(low < high, "TCP must increase with gEUD");
        // A steeper slope (larger γ50) sharpens the response around TCD50.
        let steep = tcp_logistic(80.0, tcd50, 4.0);
        assert!(
            steep > high,
            "steeper γ50 must raise TCP further above TCD50"
        );
    }

    #[test]
    fn ntcp_lkb_matches_the_normal_cdf() {
        // At gEUD = TD50, t = 0 ⇒ Φ(0) = 0.5.
        let (td50, m) = (50.0, 0.2);
        assert_relative_eq!(ntcp_lkb(td50, td50, m), 0.5, epsilon = 1e-12);
        // gEUD = 60 ⇒ t = (60−50)/(0.2·50) = 1 ⇒ Φ(1) = 0.841344746 (published).
        assert_relative_eq!(ntcp_lkb(60.0, td50, m), 0.841_344_746, max_relative = 1e-6);
        // gEUD = 40 ⇒ t = −1 ⇒ Φ(−1) = 1 − Φ(1) = 0.158655254.
        assert_relative_eq!(ntcp_lkb(40.0, td50, m), 0.158_655_254, max_relative = 1e-5);
        // Bounded and monotone.
        assert!(ntcp_lkb(30.0, td50, m) < ntcp_lkb(70.0, td50, m));
        assert!((0.0..=1.0).contains(&ntcp_lkb(200.0, td50, m)));
    }

    #[test]
    fn radiobiology_is_generic_over_scalar_f32() {
        assert_relative_eq!(generalized_eud(&[2.0_f32, 4.0], 1.0), 3.0, epsilon = 1e-6);
        assert_relative_eq!(tcp_logistic(50.0_f32, 50.0, 2.0), 0.5, epsilon = 1e-6);
        assert_relative_eq!(ntcp_lkb(50.0_f32, 50.0, 0.2), 0.5, epsilon = 1e-6);
    }
}
