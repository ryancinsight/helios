//! Kernel superposition — stage 2 of the collapsed-cone / convolution dose model.
//!
//! Stage 1 ([`deposit_ray_terma`](crate::deposit_ray_terma),
//! [`primary_fluence_parallel_x`](crate::primary_fluence_parallel_x)) yields the
//! **terma** — total energy released per voxel by the primary beam. Stage 2 spreads
//! that released energy to surrounding voxels with a **dose-deposition kernel**,
//! turning terma into dose: this is where lateral penumbra and depth build-up come
//! from (the primary-only terma has neither).
//!
//! The kernel here is **separable** (`K = kₓ ⊗ k_y ⊗ k_z`), so the 3-D convolution
//! factors into three cheap axis passes (`O(N·taps)` each instead of `O(N·taps³)`).
//! A separable symmetric kernel models isotropic-ish scatter with per-axis ranges;
//! a full anisotropic (forward-peaked, poly-energetic) collapsed-cone kernel is a
//! later increment. Each axis kernel is **centred** (index `len/2` is offset 0) and
//! normalized to `Σ = 1`, so energy is conserved in the interior; a `[1]` kernel on
//! every axis is the identity (dose = terma), the differential oracle against the
//! primary-only reference.

use aequitas::systems::si::quantities::{Dimensionless, Length};
use helios_core::constants::CM_PER_M;
use helios_domain::{Volume, VoxelGrid};
use helios_math::{NumericElement, Scalar};

#[inline]
fn normalize_weights_by_sum<T: Scalar>(weights: &mut [T], sum: T) {
    let zero = <T as NumericElement>::ZERO;
    if sum > zero {
        let inv_sum = sum.recip();
        for weight in weights {
            *weight *= inv_sum;
        }
    }
}

#[inline]
fn normalized_kernel_from_taps<T, F>(taps: usize, mut tap_weight: F) -> Vec<T>
where
    T: Scalar,
    F: FnMut(usize) -> T,
{
    let zero = <T as NumericElement>::ZERO;
    let mut kernel = Vec::with_capacity(taps);
    let mut sum = zero;
    for tap in 0..taps {
        let weight = tap_weight(tap);
        kernel.push(weight);
        sum += weight;
    }
    normalize_weights_by_sum(&mut kernel, sum);
    kernel
}

/// Symmetric normalized deposition kernel `k[d] ∝ exp(−|offset|·voxel_spacing / range)`
/// over offsets `[−radius, radius]` (length `2·radius + 1`), normalized so `Σ = 1`.
///
/// `range` is the characteristic scatter/transport range; `voxel_spacing` is the
/// voxel spacing along the axis. Both values retain their Aequitas units through
/// this API and are converted to centimetres only at the exponential formula
/// boundary. `radius = 0` returns `[1]` (the identity / no-spread kernel).
#[must_use]
pub fn symmetric_deposition_kernel<T: Scalar>(
    range: Length<T>,
    voxel_spacing: Length<T>,
    radius: usize,
) -> Vec<T> {
    let taps = 2 * radius + 1;
    let range_cm = range.into_base() * T::from_f64(CM_PER_M);
    let voxel_cm = voxel_spacing.into_base() * T::from_f64(CM_PER_M);
    let inv_range = range_cm.recip();
    normalized_kernel_from_taps(taps, |t| {
        let offset = (t as f64 - radius as f64).abs();
        let distance = T::from_f64(offset) * voxel_cm;
        (-(distance * inv_range)).exp()
    })
}

/// Convolve `vol` with a centred 1-D `kernel` along `axis` (0 = x, 1 = y, 2 = z).
///
/// Offset-0 is kernel index `len/2`; taps whose source voxel falls outside the grid
/// are dropped (energy leaving the boundary is not wrapped), so interior voxels are
/// exact while the boundary layer loses the truncated tail.
///
/// Iterates the volume's zero-copy [`as_slice`](Volume::as_slice) view with a
/// compile-time axis and precomputed stride. Each voxel derives its contiguous
/// valid tap interval once, eliminating per-tap boundary branches while retaining
/// the original tap summation order. Results are therefore bitwise-identical to
/// the bounds-checked reference form.
fn convolve_axis<T: Scalar, const AXIS: usize>(vol: &Volume<T>, kernel: &[T]) -> Volume<T> {
    convolve_axis_at_const::<T, AXIS>(vol, kernel, kernel.len() / 2)
}

/// [`convolve_axis`] with an explicit zero-offset index `center` — the general
/// form serving **asymmetric** (forward-peaked) kernels, where offset 0 is not
/// the midpoint. `center = len/2` recovers the centred behaviour exactly.
fn convolve_axis_at<T: Scalar>(
    vol: &Volume<T>,
    kernel: &[T],
    center: usize,
    axis: usize,
) -> Volume<T> {
    match axis {
        0 => convolve_axis_at_const::<T, 0>(vol, kernel, center),
        1 => convolve_axis_at_const::<T, 1>(vol, kernel, center),
        2 => convolve_axis_at_const::<T, 2>(vol, kernel, center),
        _ => panic!("axis must be 0, 1, or 2; got {axis}"),
    }
}

fn convolve_axis_at_const<T: Scalar, const AXIS: usize>(
    vol: &Volume<T>,
    kernel: &[T],
    center: usize,
) -> Volume<T> {
    const { assert!(AXIS < 3, "axis must be 0, 1, or 2") };
    let grid: VoxelGrid<T> = *vol.grid();
    let [nx, ny, nz] = grid.dims();
    let extent = [nx, ny, nz][AXIS];
    // C-contiguous (i, j, k) strides — the Volume layout contract.
    let axis_stride = [ny * nz, nz, 1][AXIS] as isize;

    let src = vol.as_slice();
    let mut out = vec![<T as NumericElement>::ZERO; src.len()];
    let mut base = 0usize;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let pos = [i, j, k][AXIS];
                let centered_pos = pos + center;
                let first_tap = centered_pos.saturating_add(1).saturating_sub(extent);
                let past_last_tap = centered_pos.saturating_add(1).min(kernel.len());
                let mut acc = <T as NumericElement>::ZERO;
                for (relative_tap, &weight) in kernel[first_tap..past_last_tap].iter().enumerate() {
                    // True convolution: dose(pos) gathers src(pos − offset), so a
                    // downstream-weighted (offset > 0) tap carries energy FROM the
                    // upstream source TO pos. (Correlation — pos + offset — would
                    // invert asymmetric kernels; symmetric ones cannot tell.)
                    let source_pos = centered_pos - (first_tap + relative_tap);
                    let offset = base as isize + (source_pos as isize - pos as isize) * axis_stride;
                    acc += src[offset as usize] * weight;
                }
                out[base] = acc;
                base += 1;
            }
        }
    }
    Volume::from_shape_vec(grid, out).expect("output length equals input voxel count")
}

/// Dose by separable 3-D convolution-superposition of a `terma` volume with
/// centred per-axis deposition kernels `kx`, `ky`, `kz`.
///
/// Applies the three axis convolutions in turn. With `[1]` kernels this is the
/// identity (`dose = terma`); with normalized spread kernels it reproduces lateral
/// penumbra (energy from a beamlet reaches neighbouring voxels) and, along the beam,
/// depth build-up — and conserves energy in the interior. Linear in `terma`.
#[must_use]
pub fn scatter_superposition<T: Scalar>(
    terma: &Volume<T>,
    kx: &[T],
    ky: &[T],
    kz: &[T],
) -> Volume<T> {
    let after_x = convolve_axis::<T, 0>(terma, kx);
    let after_y = convolve_axis::<T, 1>(&after_x, ky);
    convolve_axis::<T, 2>(&after_y, kz)
}

/// Forward-peaked (anisotropic) deposition kernel along the beam axis:
/// `k[d] ∝ exp(−|offset|·voxel_spacing / range)` with **different ranges upstream
/// vs downstream** — `range_down` (beam direction, secondary-electron forward
/// transport) and `range_up` (backscatter, physically much shorter). Offsets
/// span `[−radius_up, +radius_down]`; the returned `usize` is the zero-offset
/// index (`radius_up`). Normalized so `Σ = 1` (energy-conserving in the
/// interior). Equal ranges and radii reduce to
/// [`symmetric_deposition_kernel`] exactly (the differential oracle).
#[must_use]
pub fn forward_peaked_kernel<T: Scalar>(
    range_up: Length<T>,
    range_down: Length<T>,
    voxel_spacing: Length<T>,
    radius_up: usize,
    radius_down: usize,
) -> (Vec<T>, usize) {
    let taps = radius_up + radius_down + 1;
    let range_up_cm = range_up.into_base() * T::from_f64(CM_PER_M);
    let range_down_cm = range_down.into_base() * T::from_f64(CM_PER_M);
    let voxel_cm = voxel_spacing.into_base() * T::from_f64(CM_PER_M);
    let (inv_up, inv_down) = (range_up_cm.recip(), range_down_cm.recip());
    let kernel = normalized_kernel_from_taps(taps, |t| {
        let offset = t as f64 - radius_up as f64; // <0 upstream, >0 downstream
        let distance = T::from_f64(offset.abs()) * voxel_cm;
        let inv_range = if offset < 0.0 { inv_up } else { inv_down };
        (-(distance * inv_range)).exp()
    });
    (kernel, radius_up)
}

/// One spectral component of a poly-energetic beam for
/// [`poly_forward_peaked_kernel`]: its forward-peaked upstream/downstream
/// transport ranges and its relative energy-fluence `weight`.
///
/// Higher-energy components carry farther downstream (larger `range_down`),
/// so a spectrum weighted toward high energy is more forward-peaked — the
/// beam-hardening signature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralComponent<T: Scalar> {
    /// Upstream (backscatter) transport range.
    pub range_up: Length<T>,
    /// Downstream (forward) transport range.
    pub range_down: Length<T>,
    /// Relative energy-fluence weight (need not be pre-normalized).
    pub weight: Dimensionless<T>,
}

/// Poly-energetic forward-peaked deposition kernel: the energy-fluence-weighted
/// convex combination of the monoenergetic [`forward_peaked_kernel`]s of each
/// `components` entry (all sharing `radius_up`/`radius_down`, so they superpose
/// tap-for-tap). Models beam hardening — a real MV beam is a spectrum, and its
/// harder components reach farther downstream.
///
/// Because each monoenergetic kernel is already `Σ = 1`, the weighted sum sums to
/// the total weight and is renormalized to `Σ = 1` here, so `weight`s need not be
/// pre-normalized (the result is scale-invariant in the weights). A single
/// positive-weight component reduces **exactly** to [`forward_peaked_kernel`]
/// (the differential oracle). With no positive-weight component the kernel is the
/// centred delta (identity — no spread). Returns the kernel and its zero-offset
/// index (`radius_up`).
#[must_use]
pub fn poly_forward_peaked_kernel<T: Scalar>(
    components: &[SpectralComponent<T>],
    voxel_spacing: Length<T>,
    radius_up: usize,
    radius_down: usize,
) -> (Vec<T>, usize) {
    let zero = <T as NumericElement>::ZERO;
    let taps = radius_up + radius_down + 1;
    let mut acc = vec![zero; taps];
    let mut total_weight = zero;
    for component in components {
        let weight = component.weight.into_base();
        if weight <= zero {
            continue; // non-positive weight contributes nothing.
        }
        let (mono, _) = forward_peaked_kernel(
            component.range_up,
            component.range_down,
            voxel_spacing,
            radius_up,
            radius_down,
        );
        for (a, &m) in acc.iter_mut().zip(&mono) {
            *a += m * weight;
        }
        total_weight += weight;
    }
    if total_weight > zero {
        normalize_weights_by_sum(&mut acc, total_weight);
    } else {
        acc[radius_up] = <T as NumericElement>::ONE; // degenerate ⇒ identity.
    }
    (acc, radius_up)
}

/// Dose by **beam-aligned anisotropic** separable superposition: the
/// forward-peaked `(beam_kernel, beam_center)` (from [`forward_peaked_kernel`])
/// applies along `beam_axis` (0 = x, 1 = y, 2 = z) and the symmetric `lateral`
/// kernel along the two remaining axes.
///
/// This is the collapsed-cone anisotropy for an axis-aligned beam: more energy
/// carried downstream than upstream (build-up/downstream tail), symmetric
/// penumbra laterally. With equal up/down ranges it reduces **exactly** to
/// [`scatter_superposition`] — the differential oracle. Rotated (per-gantry-
/// angle) cone axes remain future work; the helical geometry applies this in
/// the beam's eye view.
#[must_use]
pub fn anisotropic_scatter_superposition<T: Scalar>(
    terma: &Volume<T>,
    beam_axis: usize,
    beam_kernel: &[T],
    beam_center: usize,
    lateral: &[T],
) -> Volume<T> {
    let mut vol = convolve_axis_at(terma, beam_kernel, beam_center, beam_axis);
    for axis in 0..3 {
        if axis != beam_axis {
            vol = convolve_axis_at(&vol, lateral, lateral.len() / 2, axis);
        }
    }
    vol
}

#[cfg(test)]
mod tests;
