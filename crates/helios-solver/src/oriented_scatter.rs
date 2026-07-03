//! Oriented (beam-frame) collapsed-cone scatter along an arbitrary direction.
//!
//! [`scatter`](crate::scatter) convolves along grid axes — physically correct
//! only when the beam is grid-aligned. Helical TomoTherapy rotates the gantry,
//! so each frame's beam travels an **oblique, in-plane** direction. This module
//! convolves the forward-peaked deposition kernel along an arbitrary **unit**
//! direction by trilinearly resampling the terma field, giving a beam-frame
//! collapsed cone that follows the gantry.
//!
//! The primitive is [`directional_convolve`] (an oriented 1-D convolution along
//! any unit vector); [`oriented_forward_scatter`] composes it into a full
//! collapsed cone (forward-peaked along the beam, symmetric across the two
//! lateral directions). When the beam is a grid axis and the sample step equals
//! that axis's voxel pitch, sample points land on nodes and the result matches
//! [`anisotropic_scatter_superposition`](crate::anisotropic_scatter_superposition)
//! up to trilinear-at-node exactness (the differential oracle).

use helios_domain::{Volume, VoxelGrid};
use helios_math::{GeometryScalar, NumericElement, Point3, Vector3};

/// Convolve `vol` with the 1-D `kernel` (zero-offset index `center`), sampled
/// along the unit `direction` at `sample_step_mm` spacing, by trilinear gather:
///
/// ```text
/// out(p) = Σ_t kernel[t] · vol.sample_world(p − (t − center)·step·direction)
/// ```
///
/// Samples that fall outside the grid contribute `0` (boundary truncation,
/// matching the axis convolution). This is **convolution**, not correlation: a
/// downstream-weighted tap (signed offset `t − center > 0`) gathers from the
/// upstream source, carrying energy forward along `direction` — the same sign
/// convention as [`scatter`](crate::scatter).
///
/// `direction` is assumed unit-length; pass the beam's forward unit vector.
#[must_use]
pub fn directional_convolve<T: GeometryScalar>(
    vol: &Volume<T>,
    kernel: &[T],
    center: usize,
    direction: Vector3<T>,
    sample_step_mm: T,
) -> Volume<T> {
    let grid: VoxelGrid<T> = *vol.grid();
    let [nx, ny, nz] = grid.dims();
    let zero = <T as NumericElement>::ZERO;
    let center_f = <T as GeometryScalar>::from_f64(center as f64);

    let mut out = vec![zero; nx * ny * nz];
    let mut idx = 0usize;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let p = grid.voxel_center(i, j, k);
                let mut acc = zero;
                for (t, &weight) in kernel.iter().enumerate() {
                    let signed = <T as GeometryScalar>::from_f64(t as f64) - center_f;
                    let shift = signed * sample_step_mm;
                    let sample_point = Point3::new(
                        p.x - direction.x * shift,
                        p.y - direction.y * shift,
                        p.z - direction.z * shift,
                    );
                    if let Some(v) = vol.sample_world(sample_point) {
                        acc += v * weight;
                    }
                }
                out[idx] = acc;
                idx += 1;
            }
        }
    }
    Volume::from_shape_vec(grid, out).expect("output length equals input voxel count")
}

/// Orthonormal lateral basis `(u, v)` spanning the plane perpendicular to the
/// unit `beam` direction.
///
/// Seeds Gram–Schmidt from the world axis least parallel to `beam` (so the
/// projection is well-conditioned), then `u = normalize(seed − (seed·beam)beam)`
/// and `v = beam × u`. For an in-plane beam (`beam.z = 0`) this yields
/// `u = +z` and `v` the in-plane perpendicular — the natural TomoTherapy frame.
fn lateral_basis<T: GeometryScalar>(beam: Vector3<T>) -> (Vector3<T>, Vector3<T>) {
    let zero = <T as NumericElement>::ZERO;
    let one = <T as NumericElement>::ONE;
    // |beam·ẑ| < 0.9 ⇒ ẑ is a safe seed; otherwise the beam is near-vertical, use x̂.
    let seed = if beam.z.to_f64().abs() < 0.9 {
        Vector3::new(zero, zero, one)
    } else {
        Vector3::new(one, zero, zero)
    };
    let u = (seed - beam * seed.dot(beam)).normalize();
    let v = beam.cross(u);
    (u, v)
}

/// Beam-frame collapsed-cone dose from a `terma` field: the forward-peaked
/// `(beam_kernel, beam_center)` convolved along the (arbitrary) unit `beam_dir`,
/// then the symmetric `lateral` kernel convolved across both directions
/// perpendicular to the beam. `sample_step_mm` is the physical sampling pitch
/// along every direction (the `voxel_cm` the kernels were built with, ×10).
///
/// This is the gantry-following anisotropy: more energy carried downstream than
/// upstream **along the actual beam**, symmetric penumbra laterally, at any
/// gantry angle. With `beam_dir` a grid axis and `sample_step_mm` that axis's
/// pitch it reduces to
/// [`anisotropic_scatter_superposition`](crate::anisotropic_scatter_superposition)
/// (up to trilinear-at-node exactness).
#[must_use]
pub fn oriented_forward_scatter<T: GeometryScalar>(
    terma: &Volume<T>,
    beam_dir: Vector3<T>,
    beam_kernel: &[T],
    beam_center: usize,
    lateral: &[T],
    sample_step_mm: T,
) -> Volume<T> {
    let beam = beam_dir.normalize();
    let (u, v) = lateral_basis(beam);
    let lateral_center = lateral.len() / 2;

    let after_beam = directional_convolve(terma, beam_kernel, beam_center, beam, sample_step_mm);
    let after_u = directional_convolve(&after_beam, lateral, lateral_center, u, sample_step_mm);
    directional_convolve(&after_u, lateral, lateral_center, v, sample_step_mm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        anisotropic_scatter_superposition, forward_peaked_kernel, symmetric_deposition_kernel,
    };
    use approx::assert_relative_eq;
    use helios_math::Point3 as P3;

    fn grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], P3::new(0.0, 0.0, 0.0)).expect("grid")
    }

    fn point_terma() -> Volume<f64> {
        Volume::from_shape_fn(grid(), |idx| if idx == [4, 4, 4] { 1.0 } else { 0.0 })
    }

    #[test]
    fn axis_aligned_reduces_to_the_separable_anisotropic_scatter() {
        // beam = +x, sample step = x pitch (2 mm) ⇒ every sample lands on a node ⇒
        // trilinear is exact ⇒ the oriented result must match the grid-axis
        // separable form (the differential oracle) to trilinear precision.
        let (fp, centre) = forward_peaked_kernel(0.15_f64, 0.6, 0.2, 1, 3);
        let lat = symmetric_deposition_kernel(0.4_f64, 0.2, 1);
        let terma = point_terma();

        let separable = anisotropic_scatter_superposition(&terma, 0, &fp, centre, &lat);
        let oriented =
            oriented_forward_scatter(&terma, Vector3::new(1.0, 0.0, 0.0), &fp, centre, &lat, 2.0);

        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    assert_relative_eq!(
                        oriented.get(i, j, k).unwrap(),
                        separable.get(i, j, k).unwrap(),
                        epsilon = 1e-10
                    );
                }
            }
        }
    }

    #[test]
    fn oblique_beam_deposits_more_energy_downstream_than_upstream() {
        // 45° in-plane beam. From the point source, sample the dose one step
        // downstream (+beam) vs upstream (−beam) in world space: forward-peaking
        // must put strictly more dose downstream, and the two lateral flanks
        // (±perp) must stay symmetric.
        let (fp, centre) = forward_peaked_kernel(0.1_f64, 1.0, 0.2, 1, 3);
        let lat = symmetric_deposition_kernel(0.3_f64, 0.2, 1);
        let inv = 1.0 / 2.0_f64.sqrt();
        let beam = Vector3::new(inv, inv, 0.0);
        let dose = oriented_forward_scatter(&point_terma(), beam, &fp, centre, &lat, 2.0);

        let src = grid().voxel_center(4, 4, 4);
        let d = 4.0; // mm along the beam
        let down = dose
            .sample_world(P3::new(src.x + beam.x * d, src.y + beam.y * d, src.z))
            .unwrap();
        let up = dose
            .sample_world(P3::new(src.x - beam.x * d, src.y - beam.y * d, src.z))
            .unwrap();
        assert!(down > up, "downstream {down} must exceed upstream {up}");

        // Lateral flanks along the in-plane perpendicular (−y, +x)/√2 are symmetric.
        let perp = Vector3::new(-inv, inv, 0.0);
        let left = dose
            .sample_world(P3::new(src.x + perp.x * d, src.y + perp.y * d, src.z))
            .unwrap();
        let right = dose
            .sample_world(P3::new(src.x - perp.x * d, src.y - perp.y * d, src.z))
            .unwrap();
        assert_relative_eq!(left, right, epsilon = 1e-12);
    }

    #[test]
    fn oblique_scatter_approximately_conserves_interior_energy() {
        // Σ=1 kernels + trilinear gather conserve interior energy up to boundary
        // truncation and interpolation diffusion. The source sits at the grid
        // centre, far from every face; assert conservation within 2 %.
        let (fp, centre) = forward_peaked_kernel(0.2_f64, 0.5, 0.2, 1, 2);
        let lat = symmetric_deposition_kernel(0.3_f64, 0.2, 1);
        let inv = 1.0 / 2.0_f64.sqrt();
        let dose = oriented_forward_scatter(
            &point_terma(),
            Vector3::new(inv, inv, 0.0),
            &fp,
            centre,
            &lat,
            2.0,
        );
        assert_relative_eq!(dose.sum(), 1.0, max_relative = 0.02);
    }

    #[test]
    fn oriented_scatter_is_generic_over_scalar_f32() {
        let g = VoxelGrid::<f32>::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], P3::new(0.0, 0.0, 0.0))
            .unwrap();
        let terma = Volume::from_shape_fn(g, |idx| if idx == [4, 4, 4] { 1.0_f32 } else { 0.0 });
        let (fp, centre) = forward_peaked_kernel(0.1_f32, 1.0, 0.2, 1, 3);
        let lat = symmetric_deposition_kernel(0.3_f32, 0.2, 1);
        let inv = 1.0_f32 / 2.0_f32.sqrt();
        let beam = Vector3::new(inv, inv, 0.0);
        let dose = oriented_forward_scatter(&terma, beam, &fp, centre, &lat, 2.0);
        let src = g.voxel_center(4, 4, 4);
        let down = dose
            .sample_world(P3::new(src.x + beam.x * 4.0, src.y + beam.y * 4.0, src.z))
            .unwrap();
        let up = dose
            .sample_world(P3::new(src.x - beam.x * 4.0, src.y - beam.y * 4.0, src.z))
            .unwrap();
        assert!(down > up);
    }
}
