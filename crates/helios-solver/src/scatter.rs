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

use helios_domain::{Volume, VoxelGrid};
use helios_math::{NumericElement, Scalar};

/// Symmetric normalized deposition kernel `k[d] ∝ exp(−|offset|·voxel_cm / range_cm)`
/// over offsets `[−radius, radius]` (length `2·radius + 1`), normalized so `Σ = 1`.
///
/// `range_cm` is the characteristic scatter/transport range; `voxel_cm` the voxel
/// spacing along the axis. `radius = 0` returns `[1]` (the identity / no-spread
/// kernel).
#[must_use]
pub fn symmetric_deposition_kernel<T: Scalar>(range_cm: T, voxel_cm: T, radius: usize) -> Vec<T> {
    let zero = <T as NumericElement>::ZERO;
    let taps = 2 * radius + 1;
    let inv_range = range_cm.recip();
    let mut kernel = Vec::with_capacity(taps);
    let mut sum = zero;
    for t in 0..taps {
        let offset = (t as f64 - radius as f64).abs();
        let distance = T::from_f64(offset) * voxel_cm;
        let weight = (-(distance * inv_range)).exp();
        kernel.push(weight);
        sum += weight;
    }
    if sum > zero {
        let inv_sum = sum.recip();
        for w in &mut kernel {
            *w *= inv_sum;
        }
    }
    kernel
}

/// Convolve `vol` with a centred 1-D `kernel` along `axis` (0 = x, 1 = y, 2 = z).
///
/// Offset-0 is kernel index `len/2`; taps whose source voxel falls outside the grid
/// are dropped (energy leaving the boundary is not wrapped), so interior voxels are
/// exact while the boundary layer loses the truncated tail.
///
/// Iterates the volume's zero-copy [`as_slice`](Volume::as_slice) view with a
/// precomputed axis stride: per tap this is one strided slice read instead of a
/// three-axis bounds check + flat-index recomputation through `get`. The tap
/// summation order is unchanged, so results are bitwise-identical to the previous
/// per-voxel form (the differential guarantee for this optimization).
fn convolve_axis<T: Scalar>(vol: &Volume<T>, kernel: &[T], axis: usize) -> Volume<T> {
    convolve_axis_at(vol, kernel, kernel.len() / 2, axis)
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
    let grid: VoxelGrid<T> = *vol.grid();
    let [nx, ny, nz] = grid.dims();
    let center = center as isize;
    let extent = [nx, ny, nz][axis] as isize;
    // C-contiguous (i, j, k) strides — the Volume layout contract.
    let axis_stride = [ny * nz, nz, 1][axis] as isize;

    let src = vol.as_slice();
    let mut out = vec![<T as NumericElement>::ZERO; src.len()];
    let mut base = 0usize;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let pos = [i, j, k][axis] as isize;
                let mut acc = <T as NumericElement>::ZERO;
                for (t, &weight) in kernel.iter().enumerate() {
                    // True convolution: dose(pos) gathers src(pos − offset), so a
                    // downstream-weighted (offset > 0) tap carries energy FROM the
                    // upstream source TO pos. (Correlation — pos + offset — would
                    // invert asymmetric kernels; symmetric ones cannot tell.)
                    let s = pos - (t as isize - center);
                    if s < 0 || s >= extent {
                        continue;
                    }
                    let offset = base as isize + (s - pos) * axis_stride;
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
    let after_x = convolve_axis(terma, kx, 0);
    let after_y = convolve_axis(&after_x, ky, 1);
    convolve_axis(&after_y, kz, 2)
}

/// Forward-peaked (anisotropic) deposition kernel along the beam axis:
/// `k[d] ∝ exp(−|offset|·voxel_cm / range)` with **different ranges upstream vs
/// downstream** — `range_down_cm` (beam direction, secondary-electron forward
/// transport) and `range_up_cm` (backscatter, physically much shorter). Offsets
/// span `[−radius_up, +radius_down]`; the returned `usize` is the zero-offset
/// index (`radius_up`). Normalized so `Σ = 1` (energy-conserving in the
/// interior). Equal ranges and radii reduce to
/// [`symmetric_deposition_kernel`] exactly (the differential oracle).
#[must_use]
pub fn forward_peaked_kernel<T: Scalar>(
    range_up_cm: T,
    range_down_cm: T,
    voxel_cm: T,
    radius_up: usize,
    radius_down: usize,
) -> (Vec<T>, usize) {
    let zero = <T as NumericElement>::ZERO;
    let taps = radius_up + radius_down + 1;
    let (inv_up, inv_down) = (range_up_cm.recip(), range_down_cm.recip());
    let mut kernel = Vec::with_capacity(taps);
    let mut sum = zero;
    for t in 0..taps {
        let offset = t as f64 - radius_up as f64; // <0 upstream, >0 downstream
        let distance = T::from_f64(offset.abs()) * voxel_cm;
        let inv_range = if offset < 0.0 { inv_up } else { inv_down };
        let weight = (-(distance * inv_range)).exp();
        kernel.push(weight);
        sum += weight;
    }
    if sum > zero {
        let inv_sum = sum.recip();
        for w in &mut kernel {
            *w *= inv_sum;
        }
    }
    (kernel, radius_up)
}

/// One spectral component of a poly-energetic beam for
/// [`poly_forward_peaked_kernel`]: its forward-peaked upstream/downstream
/// transport ranges (cm) and its relative energy-fluence `weight`.
///
/// Higher-energy components carry farther downstream (larger `range_down_cm`),
/// so a spectrum weighted toward high energy is more forward-peaked — the
/// beam-hardening signature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralComponent<T: Scalar> {
    /// Upstream (backscatter) transport range (cm).
    pub range_up_cm: T,
    /// Downstream (forward) transport range (cm).
    pub range_down_cm: T,
    /// Relative energy-fluence weight (need not be pre-normalized).
    pub weight: T,
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
    voxel_cm: T,
    radius_up: usize,
    radius_down: usize,
) -> (Vec<T>, usize) {
    let zero = <T as NumericElement>::ZERO;
    let taps = radius_up + radius_down + 1;
    let mut acc = vec![zero; taps];
    let mut total_weight = zero;
    for component in components {
        if component.weight <= zero {
            continue; // non-positive weight contributes nothing.
        }
        let (mono, _) = forward_peaked_kernel(
            component.range_up_cm,
            component.range_down_cm,
            voxel_cm,
            radius_up,
            radius_down,
        );
        for (a, &m) in acc.iter_mut().zip(&mono) {
            *a += m * component.weight;
        }
        total_weight += component.weight;
    }
    if total_weight > zero {
        let inv = total_weight.recip();
        for a in &mut acc {
            *a *= inv;
        }
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
            vol = convolve_axis(&vol, lateral, axis);
        }
    }
    vol
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use helios_math::Point3;

    fn grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([7, 7, 7], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid")
    }

    // Terma concentrated in the single centre voxel (3,3,3).
    fn point_terma() -> Volume<f64> {
        Volume::from_shape_fn(grid(), |idx| if idx == [3, 3, 3] { 1.0 } else { 0.0 })
    }

    #[test]
    fn forward_peaked_reduces_to_symmetric_for_equal_ranges() {
        // Equal up/down ranges and radii → identical taps and centre to the
        // symmetric kernel (differential oracle for the asymmetric constructor).
        let sym = symmetric_deposition_kernel(0.5_f64, 0.2, 2);
        let (fp, centre) = forward_peaked_kernel(0.5_f64, 0.5, 0.2, 2, 2);
        assert_eq!(centre, 2);
        assert_eq!(fp.len(), sym.len());
        for (a, b) in fp.iter().zip(&sym) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
        // And the anisotropic superposition equals the symmetric one exactly.
        let iso = scatter_superposition(&point_terma(), &sym, &sym, &sym);
        let aniso = anisotropic_scatter_superposition(&point_terma(), 0, &fp, centre, &sym);
        for i in 0..7 {
            for j in 0..7 {
                for k in 0..7 {
                    assert_relative_eq!(
                        aniso.get(i, j, k).unwrap(),
                        iso.get(i, j, k).unwrap(),
                        epsilon = 1e-15
                    );
                }
            }
        }
    }

    #[test]
    fn forward_peaking_puts_more_energy_downstream_than_upstream() {
        // Long downstream range, short upstream: the voxel one step downstream
        // (+x) of the point source must receive strictly more dose than one step
        // upstream (−x) — the defining collapsed-cone anisotropy. Laterally the
        // spread stays symmetric.
        let (fp, centre) = forward_peaked_kernel(0.1_f64, 1.0, 0.2, 1, 3);
        let lat = symmetric_deposition_kernel(0.3_f64, 0.2, 1);
        let dose = anisotropic_scatter_superposition(&point_terma(), 0, &fp, centre, &lat);
        let down = dose.get(4, 3, 3).unwrap(); // +x of the source at (3,3,3)
        let up = dose.get(2, 3, 3).unwrap(); // −x
        assert!(down > up, "downstream {down} must exceed upstream {up}");
        // Lateral symmetry is preserved.
        assert_relative_eq!(
            dose.get(3, 2, 3).unwrap(),
            dose.get(3, 4, 3).unwrap(),
            epsilon = 1e-14
        );
    }

    #[test]
    fn anisotropic_kernel_conserves_interior_point_energy() {
        // Source far enough from every boundary that no tail truncates: the total
        // dose equals the total terma (Σ=1 normalization on every axis).
        let (fp, centre) = forward_peaked_kernel(0.2_f64, 0.6, 0.2, 1, 2);
        let lat = symmetric_deposition_kernel(0.4_f64, 0.2, 1);
        let dose = anisotropic_scatter_superposition(&point_terma(), 0, &fp, centre, &lat);
        assert_relative_eq!(dose.sum(), 1.0, epsilon = 1e-13);
    }

    fn spectral(range_up_cm: f64, range_down_cm: f64, weight: f64) -> SpectralComponent<f64> {
        SpectralComponent {
            range_up_cm,
            range_down_cm,
            weight,
        }
    }

    #[test]
    fn single_component_poly_reduces_to_the_monoenergetic_kernel() {
        // A one-component spectrum (any positive weight) is the monoenergetic
        // kernel exactly — the weight cancels in the Σ=1 renormalization.
        let (mono, mc) = forward_peaked_kernel(0.15_f64, 0.6, 0.2, 1, 3);
        let (poly, pc) = poly_forward_peaked_kernel(&[spectral(0.15, 0.6, 3.7)], 0.2, 1, 3);
        assert_eq!(mc, pc);
        assert_eq!(mono.len(), poly.len());
        for (a, b) in poly.iter().zip(&mono) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
    }

    #[test]
    fn poly_kernel_is_invariant_to_weight_scaling_and_sums_to_one() {
        // Scaling every weight by a constant leaves the convex combination
        // unchanged, and the kernel is normalized.
        let comps = [spectral(0.1, 0.4, 1.0), spectral(0.2, 1.2, 2.0)];
        let scaled = [spectral(0.1, 0.4, 10.0), spectral(0.2, 1.2, 20.0)];
        let (ka, _) = poly_forward_peaked_kernel(&comps, 0.2, 1, 3);
        let (kb, _) = poly_forward_peaked_kernel(&scaled, 0.2, 1, 3);
        for (a, b) in ka.iter().zip(&kb) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
        assert_relative_eq!(ka.iter().sum::<f64>(), 1.0, epsilon = 1e-14);
    }

    #[test]
    fn harder_spectrum_is_more_forward_peaked() {
        // Downstream/upstream mass ratio grows as weight shifts to the harder
        // (longer downstream range) component — the beam-hardening signature.
        let soft = spectral(0.2, 0.3, 1.0); // near-isotropic
        let hard = spectral(0.05, 1.5, 1.0); // strongly forward
        let ratio = |k: &[f64], centre: usize| -> f64 {
            let down: f64 = k[centre + 1..].iter().sum();
            let up: f64 = k[..centre].iter().sum();
            down / up
        };
        let (mostly_soft, c) = poly_forward_peaked_kernel(
            &[spectral(0.2, 0.3, 9.0), spectral(0.05, 1.5, 1.0)],
            0.2,
            2,
            2,
        );
        let (mostly_hard, _) = poly_forward_peaked_kernel(
            &[spectral(0.2, 0.3, 1.0), spectral(0.05, 1.5, 9.0)],
            0.2,
            2,
            2,
        );
        let _ = (soft, hard);
        assert!(
            ratio(&mostly_hard, c) > ratio(&mostly_soft, c),
            "harder spectrum ratio {} must exceed softer {}",
            ratio(&mostly_hard, c),
            ratio(&mostly_soft, c)
        );
    }

    #[test]
    fn empty_spectrum_is_the_identity_kernel() {
        let (k, centre) = poly_forward_peaked_kernel::<f64>(&[], 0.2, 2, 1);
        assert_eq!(centre, 2);
        assert_eq!(k, vec![0.0, 0.0, 1.0, 0.0]); // centred delta at radius_up = 2
    }

    #[test]
    fn poly_kernel_is_generic_over_scalar_f32() {
        let (mono, _) = forward_peaked_kernel(0.1_f32, 0.5, 0.2, 1, 2);
        let (poly, _) = poly_forward_peaked_kernel(
            &[SpectralComponent {
                range_up_cm: 0.1_f32,
                range_down_cm: 0.5,
                weight: 2.0,
            }],
            0.2,
            1,
            2,
        );
        for (a, b) in poly.iter().zip(&mono) {
            assert_relative_eq!(a, b, epsilon = 1e-6);
        }
    }

    #[test]
    fn anisotropic_is_generic_over_scalar_f32_and_axis_selectable() {
        let g =
            VoxelGrid::<f32>::axis_aligned([7, 7, 7], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let terma = Volume::from_shape_fn(g, |idx| if idx == [3, 3, 3] { 1.0_f32 } else { 0.0 });
        let (fp, centre) = forward_peaked_kernel(0.1_f32, 1.0, 0.2, 1, 2);
        let lat = symmetric_deposition_kernel(0.3_f32, 0.2, 1);
        // Beam along y: anisotropy shows on the j axis, not i.
        let dose = anisotropic_scatter_superposition(&terma, 1, &fp, centre, &lat);
        assert!(dose.get(3, 4, 3).unwrap() > dose.get(3, 2, 3).unwrap());
        assert_relative_eq!(
            dose.get(2, 3, 3).unwrap(),
            dose.get(4, 3, 3).unwrap(),
            epsilon = 1e-6
        );
    }

    #[test]
    fn delta_kernel_is_identity() {
        // [1] on every axis deposits energy locally → dose == terma exactly.
        let terma = Volume::from_shape_fn(grid(), |idx| (idx[0] + 2 * idx[1] + 3 * idx[2]) as f64);
        let dose = scatter_superposition(&terma, &[1.0], &[1.0], &[1.0]);
        for i in 0..7 {
            for j in 0..7 {
                for k in 0..7 {
                    assert_relative_eq!(
                        dose.get(i, j, k).unwrap(),
                        terma.get(i, j, k).unwrap(),
                        epsilon = 1e-15
                    );
                }
            }
        }
    }

    #[test]
    fn symmetric_kernel_spreads_symmetrically() {
        // Spread the point terma along x only; the two x-neighbours of the centre
        // receive equal dose (kernel symmetry), and both are < the centre.
        let kx = symmetric_deposition_kernel(0.4_f64, 0.2, 2);
        let dose = scatter_superposition(&point_terma(), &kx, &[1.0], &[1.0]);
        let left = dose.get(2, 3, 3).unwrap();
        let right = dose.get(4, 3, 3).unwrap();
        assert_relative_eq!(left, right, epsilon = 1e-14);
        assert!(left > 0.0 && left < dose.get(3, 3, 3).unwrap());
    }

    #[test]
    fn normalized_kernel_conserves_point_energy_in_interior() {
        // The centre is >= radius from every boundary, so no tail is truncated:
        // the total dose must equal the total terma (energy conservation).
        let k = symmetric_deposition_kernel(0.5_f64, 0.2, 2);
        let dose = scatter_superposition(&point_terma(), &k, &k, &k);
        assert_relative_eq!(dose.sum(), 1.0, epsilon = 1e-13);
    }

    #[test]
    fn off_axis_neighbour_receives_lateral_penumbra() {
        // A diagonal neighbour of the point source gets non-zero dose only because
        // the kernel spreads on all three axes (penumbra).
        let k = symmetric_deposition_kernel(0.5_f64, 0.2, 1);
        let dose = scatter_superposition(&point_terma(), &k, &k, &k);
        assert!(
            dose.get(4, 4, 4).unwrap() > 0.0,
            "diagonal voxel must be lit"
        );
    }

    #[test]
    fn superposition_is_linear_in_terma() {
        let k = symmetric_deposition_kernel(0.4_f64, 0.2, 2);
        let terma = point_terma();
        let scaled = Volume::from_shape_fn(*terma.grid(), |idx| {
            3.0 * terma.get(idx[0], idx[1], idx[2]).unwrap()
        });
        let d1 = scatter_superposition(&terma, &k, &k, &k);
        let d3 = scatter_superposition(&scaled, &k, &k, &k);
        assert_relative_eq!(d3.sum(), 3.0 * d1.sum(), epsilon = 1e-13);
        assert_relative_eq!(
            d3.get(3, 3, 3).unwrap(),
            3.0 * d1.get(3, 3, 3).unwrap(),
            epsilon = 1e-13
        );
    }

    #[test]
    fn symmetric_kernel_is_normalized_and_peaked_at_centre() {
        let k = symmetric_deposition_kernel(0.5_f64, 0.2, 3);
        assert_eq!(k.len(), 7);
        let sum: f64 = k.iter().sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-15);
        // Centre tap (index 3) is the maximum; symmetric about it.
        assert!(k[3] > k[2] && k[2] > k[1]);
        assert_relative_eq!(k[2], k[4], epsilon = 1e-15);
        assert_relative_eq!(k[0], k[6], epsilon = 1e-15);
    }

    #[test]
    fn scatter_is_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([5, 5, 5], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let terma = Volume::from_shape_fn(g, |idx| if idx == [2, 2, 2] { 1.0_f32 } else { 0.0 });
        let k = symmetric_deposition_kernel(0.5_f32, 0.2, 2);
        let dose = scatter_superposition(&terma, &k, &k, &k);
        // Interior point source: energy conserved.
        assert_relative_eq!(dose.sum(), 1.0_f32, epsilon = 1e-5);
    }
}
