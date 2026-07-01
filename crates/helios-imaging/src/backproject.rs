//! Shared parallel-beam back-projector.
//!
//! The geometric adjoint of the Radon transform: smear each projection reading
//! back across the reconstruction grid along its ray. Used by both filtered
//! back-projection ([`crate::filtered_back_projection`], on ramp-filtered rows
//! scaled by Δθ) and SIRT ([`crate::sirt_reconstruction`], on raw residual rows),
//! so the interpolation geometry lives in exactly one place.

use helios_domain::{Volume, VoxelGrid};
use helios_math::{GeometryScalar, NumericElement};

/// Millimetres per centimetre (detector offsets are mm, line integrals cm).
const MM_PER_CM: f64 = 10.0;

/// Back-project projection `rows` (row-major `[angle][offset]`) onto `recon`,
/// scaling every voxel by `scale`.
///
/// For each voxel the detector coordinate at each angle is computed in cm and the
/// row is sampled by linear interpolation; contributions are summed over angles.
/// The grid's axial centre is the rotation centre (matching the forward
/// projector). Requires ≥ 2 detector offsets (uniformly spaced).
pub(crate) fn back_project_rows<T: GeometryScalar>(
    angles: &[T],
    offsets: &[T],
    rows: &[T],
    recon: &VoxelGrid<T>,
    scale: T,
) -> Volume<T> {
    let zero = <T as NumericElement>::ZERO;
    let n_off = offsets.len();
    let mm_to_cm = <T as GeometryScalar>::from_f64(MM_PER_CM).recip();
    let ds_cm = (offsets[1] - offsets[0]) * mm_to_cm;
    let off0_cm = offsets[0] * mm_to_cm;
    let trig: Vec<(T, T)> = angles.iter().map(|&t| (t.cos(), t.sin())).collect();

    let grid = *recon;
    let [nx, ny, nz] = grid.dims();
    let centre = grid.voxel_center((nx - 1) / 2, (ny - 1) / 2, (nz - 1) / 2);
    let _ = nz;

    Volume::from_shape_fn(grid, |idx| {
        let world = grid.voxel_center(idx[0], idx[1], idx[2]);
        let dx = world.x - centre.x;
        let dy = world.y - centre.y;
        let mut sum = zero;
        for (a, &(cos_t, sin_t)) in trig.iter().enumerate() {
            let s_cm = (dx * cos_t + dy * sin_t) * mm_to_cm;
            let pos = (s_cm - off0_cm) * ds_cm.recip();
            let floor = pos.floor();
            let j0 = floor.to_f64() as isize;
            if j0 >= 0 && (j0 as usize) + 1 < n_off {
                let frac = pos - floor;
                let j = j0 as usize;
                let lo = rows[a * n_off + j];
                let hi = rows[a * n_off + j + 1];
                sum += lo + (hi - lo) * frac;
            }
        }
        sum * scale
    })
}
