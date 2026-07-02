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
    let grid: VoxelGrid<T> = *vol.grid();
    let [nx, ny, nz] = grid.dims();
    let center = (kernel.len() / 2) as isize;
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
                    let s = pos + (t as isize - center);
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
