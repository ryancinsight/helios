//! Parallel-beam Radon forward transform (MVCT projection sinogram).

use helios_core::HeliosError;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement, Point3, Ray, Vector3};
use helios_solver::forward_project_ray;

/// A parallel-beam sinogram: line integrals `p(θ, s)` over projection angles `θ`
/// (rad) and signed detector offsets `s` (mm from the rotation axis), stored
/// row-major `[angle][offset]`.
#[derive(Debug, Clone, PartialEq)]
pub struct Sinogram<T: GeometryScalar> {
    angles: Vec<T>,
    offsets: Vec<T>,
    data: Vec<T>,
}

impl<T: GeometryScalar> Sinogram<T> {
    /// Projection angles (rad).
    #[must_use]
    pub fn angles(&self) -> &[T] {
        &self.angles
    }

    /// Signed detector offsets (mm).
    #[must_use]
    pub fn offsets(&self) -> &[T] {
        &self.offsets
    }

    /// `(n_angles, n_offsets)`.
    #[must_use]
    pub fn dims(&self) -> (usize, usize) {
        (self.angles.len(), self.offsets.len())
    }

    /// Line integral at `(angle_index, offset_index)`.
    #[must_use]
    pub fn get(&self, angle_index: usize, offset_index: usize) -> T {
        self.data[angle_index * self.offsets.len() + offset_index]
    }

    /// Build a sinogram from explicit geometry and row-major `[angle][offset]`
    /// readings (e.g. measured detector data or a noise-perturbed projection).
    ///
    /// # Errors
    /// [`HeliosError::InvalidDomainValue`] if `data.len() != angles.len() *
    /// offsets.len()`.
    pub fn from_readings(
        angles: Vec<T>,
        offsets: Vec<T>,
        data: Vec<T>,
    ) -> Result<Self, HeliosError> {
        let expected = angles.len() * offsets.len();
        if data.len() != expected {
            return Err(HeliosError::InvalidDomainValue {
                field: "Sinogram::from_readings",
                value: data.len() as f64,
                reason: "reading count does not match angles × offsets",
            });
        }
        Ok(Self {
            angles,
            offsets,
            data,
        })
    }

    /// Apply `f` to every line-integral reading, preserving the geometry.
    ///
    /// `f` is called in row-major `[angle][offset]` order, so a stateful `f`
    /// (e.g. a seeded noise generator) is deterministic.
    #[must_use]
    pub fn map_readings<F: FnMut(T) -> T>(&self, mut f: F) -> Self {
        Self {
            angles: self.angles.clone(),
            offsets: self.offsets.clone(),
            data: self.data.iter().map(|&v| f(v)).collect(),
        }
    }
}

/// Parallel-beam Radon transform of the axial centre slice of `mu`.
///
/// For each projection angle `θ` and signed detector offset `s`, integrates `μ`
/// along the line `{ centre + s·(cosθ, sinθ) + t·(−sinθ, cosθ) }` in the axial
/// plane. `source_distance_mm` places each ray's start well outside the grid;
/// `step_mm` is the ray-march sampling step. Rays that miss record zero.
#[must_use]
pub fn parallel_beam_radon<T: GeometryScalar>(
    mu: &Volume<T>,
    angles: &[T],
    detector_offsets: &[T],
    source_distance_mm: T,
    step_mm: T,
) -> Sinogram<T> {
    let zero = <T as NumericElement>::ZERO;
    let grid = *mu.grid();
    let [nx, ny, nz] = grid.dims();
    let centre = grid.voxel_center((nx - 1) / 2, (ny - 1) / 2, (nz - 1) / 2);

    let mut data = Vec::with_capacity(angles.len() * detector_offsets.len());
    for &theta in angles {
        let (cos_t, sin_t) = (theta.cos(), theta.sin());
        // Integration direction (along the line) and the detector-offset axis.
        let dir = Vector3::new(-sin_t, cos_t, zero);
        for &s in detector_offsets {
            // Point on the line at signed distance s from the axis.
            let px = centre.x + s * cos_t;
            let py = centre.y + s * sin_t;
            let origin = Point3::new(
                px - dir.x * source_distance_mm,
                py - dir.y * source_distance_mm,
                centre.z,
            );
            let integral = Ray::try_from_direction(origin, dir)
                .and_then(|ray| forward_project_ray(mu, &ray, step_mm))
                .unwrap_or(zero);
            data.push(integral);
        }
    }

    Sinogram {
        angles: angles.to_vec(),
        offsets: detector_offsets.to_vec(),
        data,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    /// Uniform-μ disk of radius `radius_mm` centred in a 0.5 mm grid, single axial
    /// slice. The fine voxel keeps the voxelized-circle vs analytical-circle
    /// staircase error below ~1% so the sinogram tolerance can be tight.
    fn disk_phantom(mu0: f64, radius_mm: f64) -> Volume<f64> {
        let n = 161;
        let spacing = 0.5;
        let grid = VoxelGrid::axis_aligned(
            [n, n, 1],
            [spacing, spacing, spacing],
            Point3::new(0.0, 0.0, 0.0),
        )
        .unwrap();
        let centre = (n - 1) as f64 * spacing / 2.0; // world centre (mm)
        Volume::from_shape_fn(grid, move |idx| {
            let dx = idx[0] as f64 * spacing - centre;
            let dy = idx[1] as f64 * spacing - centre;
            if (dx * dx + dy * dy).sqrt() <= radius_mm {
                mu0
            } else {
                0.0
            }
        })
    }

    #[test]
    fn disk_sinogram_matches_analytical_chord() {
        // p(θ,s) = μ₀·2√(R²−s²) (mm) · 0.1 (mm→cm), independent of θ.
        let mu0 = 0.04;
        let radius = 25.0;
        let vol = disk_phantom(mu0, radius);
        let angles = [
            0.0_f64,
            std::f64::consts::FRAC_PI_4,
            std::f64::consts::FRAC_PI_2,
        ];
        let offsets = [-15.0_f64, 0.0, 10.0];
        let sino = parallel_beam_radon(&vol, &angles, &offsets, 400.0, 0.25);

        for (ai, _) in angles.iter().enumerate() {
            for (di, &s) in offsets.iter().enumerate() {
                let chord_mm = 2.0 * (radius * radius - s * s).sqrt();
                let expected = mu0 * chord_mm * 0.1;
                // ~2% tolerance: voxelized disk edge vs the analytical circle.
                assert_relative_eq!(sino.get(ai, di), expected, max_relative = 2e-2);
            }
        }
    }

    #[test]
    fn sinogram_is_angle_independent_for_a_disk() {
        let vol = disk_phantom(0.04, 25.0);
        let angles = [0.0_f64, 0.3, 1.1, 2.0];
        let offsets = [0.0_f64];
        let sino = parallel_beam_radon(&vol, &angles, &offsets, 400.0, 0.25);
        let central = sino.get(0, 0);
        for ai in 1..angles.len() {
            assert_relative_eq!(sino.get(ai, 0), central, max_relative = 2e-2);
        }
    }

    #[test]
    fn ray_outside_the_disk_reads_zero() {
        let vol = disk_phantom(0.04, 25.0);
        // Offset beyond the disk radius → the line misses the object.
        let sino = parallel_beam_radon(&vol, &[0.0_f64], &[30.0_f64], 400.0, 0.25);
        assert_relative_eq!(sino.get(0, 0), 0.0, epsilon = 1e-9);
    }

    #[test]
    fn dims_and_indexing_are_consistent() {
        let vol = disk_phantom(0.04, 20.0);
        let sino = parallel_beam_radon(&vol, &[0.0_f64, 1.0], &[-5.0_f64, 0.0, 5.0], 400.0, 0.5);
        assert_eq!(sino.dims(), (2, 3));
        assert_eq!(sino.angles().len(), 2);
        assert_eq!(sino.offsets().len(), 3);
    }

    #[test]
    fn from_readings_validates_length_and_map_preserves_geometry() {
        // Correct length constructs; wrong length errors.
        let ok = Sinogram::from_readings(
            vec![0.0_f64, 1.0],
            vec![-1.0, 1.0],
            vec![1.0, 2.0, 3.0, 4.0],
        );
        assert!(ok.is_ok());
        let bad = Sinogram::from_readings(vec![0.0_f64], vec![-1.0, 1.0], vec![1.0]);
        assert!(bad.is_err());
        // map_readings preserves geometry and applies f in order.
        let doubled = ok.unwrap().map_readings(|v| v * 2.0);
        assert_eq!(doubled.dims(), (2, 2));
        assert_eq!(doubled.get(1, 1), 8.0);
    }
}
