//! Resident GPU forward projection: batched `∫μ dl` with the volume on-device.
//!
//! The H-043b residency step: [`GpuProjector`] uploads the attenuation volume
//! **once** and then evaluates whole ray batches (a full sinogram, a delivery's
//! beamlet set) per dispatch through hephaestus's volume ray-integral kernel —
//! one rays upload and one scalar-per-ray download per batch, instead of the
//! per-op round-trips that left the elementwise path PCIe-bound (see the H-043
//! validation report). Differentially validated against the CPU
//! [`forward_project_ray`](helios_solver::forward_project_ray) reference.

use helios_core::constants::MM_PER_CM;
use helios_domain::Volume;
use helios_math::UnitQuaternion;
use hephaestus_core::{BlockWidth, ComputeDevice, HephaestusError, Result};
use hephaestus_wgpu::{ray_line_integrals_into, FieldGeometry, WgpuBuffer, WgpuDevice, RAY_STRIDE};

/// An attenuation volume resident on the GPU, ready for batched forward
/// projection.
///
/// The μ field (cm⁻¹ on a mm grid, matching the CPU projector's convention) is
/// uploaded at construction and reused across every [`Self::project_into`]
/// call — the amortization that makes the GPU path pay.
pub struct GpuProjector {
    field: WgpuBuffer<f32>,
    geometry: FieldGeometry,
}

impl GpuProjector {
    /// Upload `mu` (linear attenuation, cm⁻¹) to the device.
    ///
    /// # Errors
    /// [`HephaestusError`] on upload failure or if a grid dimension exceeds the
    /// kernel's exact-`f32` count limit (2²⁴). The current Hephaestus field
    /// kernel accepts axis-aligned geometry only, so a non-identity grid pose is
    /// rejected before upload rather than projected with discarded orientation.
    pub fn new(device: &WgpuDevice, mu: &Volume<f32>) -> Result<Self> {
        let geometry = Self::field_geometry(mu)?;
        let field = device.upload(mu.as_slice())?;
        Ok(Self { field, geometry })
    }

    /// Build the complete field metadata accepted by the current Hephaestus
    /// volume kernel.
    fn field_geometry(mu: &Volume<f32>) -> Result<FieldGeometry> {
        let grid = mu.grid();
        if grid.pose().rotation != UnitQuaternion::identity() {
            return Err(HephaestusError::DispatchFailed {
                message: "oriented voxel grids require a Hephaestus field-geometry pose kernel"
                    .into(),
            });
        }
        let [nx, ny, nz] = grid.dims();
        let dims_u32 = |n: usize, label: &str| -> Result<u32> {
            u32::try_from(n).map_err(|_| HephaestusError::DispatchFailed {
                message: format!("{label} = {n} exceeds u32 range"),
            })
        };
        let spacing = grid.spacing();
        let origin = grid.origin();
        Ok(FieldGeometry {
            dims: [
                dims_u32(nx, "dims.x")?,
                dims_u32(ny, "dims.y")?,
                dims_u32(nz, "dims.z")?,
            ],
            origin: [origin.x, origin.y, origin.z],
            spacing,
        })
    }

    /// Forward-project a packed ray batch, writing one optical depth `τ = ∫μ dl`
    /// (dimensionless: mm path lengths × cm⁻¹ μ, converted mm→cm) per ray.
    ///
    /// `rays` holds [`RAY_STRIDE`] `f32`s per ray (`origin.xyz` in mm, then a
    /// **unit** direction). Rays that miss the grid record `0` (matching the
    /// acquisition layer's miss convention). `step_mm` is the nominal midpoint
    /// sampling step, identical in meaning to the CPU projector's.
    ///
    /// # Errors
    /// [`HephaestusError`] on transfer/dispatch failure or a `rays`/`out` length
    /// mismatch.
    pub fn project_into(
        &self,
        device: &WgpuDevice,
        rays: &[f32],
        step_mm: f32,
        out: &mut [f32],
    ) -> Result<()> {
        if rays.len() != out.len() * RAY_STRIDE {
            return Err(HephaestusError::LengthMismatch {
                host_len: rays.len(),
                device_len: out.len() * RAY_STRIDE,
            });
        }
        if out.is_empty() {
            return Ok(());
        }
        let ray_buf = device.upload(rays)?;
        let out_buf = device.alloc_zeroed::<f32>(out.len())?;
        ray_line_integrals_into(
            device,
            &self.field,
            self.geometry,
            &ray_buf,
            step_mm,
            &out_buf,
            BlockWidth::DEFAULT,
        )?;
        device.download(&out_buf, out)?;
        // Kernel integrates in world (mm) units; μ is cm⁻¹ → scale mm→cm.
        let scale = (MM_PER_CM as f32).recip();
        for tau in out.iter_mut() {
            *tau *= scale;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use helios_domain::VoxelGrid;
    use helios_math::{Point3, Ray, Vector3};
    use helios_solver::forward_project_ray;

    fn mu_volume() -> Volume<f32> {
        // Heterogeneous field so the trilinear path is exercised: 17³ nodes,
        // 2 mm spacing, μ varying over all three axes.
        let grid = VoxelGrid::<f32>::axis_aligned(
            [17, 17, 17],
            [2.0, 2.0, 2.0],
            Point3::new(0.0, 0.0, 0.0),
        )
        .expect("grid");
        Volume::from_shape_fn(grid, |idx| {
            0.002 * idx[0] as f32 + 0.0015 * idx[1] as f32 + 0.001 * idx[2] as f32 + 0.02
        })
    }

    #[test]
    fn oriented_grid_is_rejected_before_gpu_upload() {
        let rotation = UnitQuaternion::try_from_rotation_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            1.0e-5,
        )
        .expect("right-handed quarter-turn basis");
        let grid = VoxelGrid::oriented(
            [2, 2, 2],
            [1.0, 1.0, 1.0],
            Point3::new(0.0, 0.0, 0.0),
            rotation,
        )
        .expect("grid");
        let mu = Volume::from_shape_fn(grid, |_| 0.02);

        let error = GpuProjector::field_geometry(&mu).expect_err("pose must be represented");
        assert!(matches!(
            error,
            HephaestusError::DispatchFailed { ref message }
                if message == "oriented voxel grids require a Hephaestus field-geometry pose kernel"
        ));
    }

    #[test]
    fn gpu_batch_matches_cpu_projector_per_ray() {
        let Ok(device) = crate::default_device() else {
            eprintln!("no GPU adapter — skipping GPU projection test");
            return;
        };
        let mu = mu_volume();
        let projector = GpuProjector::new(&device, &mu).expect("upload");

        // A fan of rays: axial rotations + one oblique + one miss.
        let mut rays = Vec::new();
        let mut dirs = Vec::new();
        for a in 0..8 {
            let theta = a as f32 * std::f32::consts::PI / 8.0;
            let (dx, dy) = (theta.cos(), theta.sin());
            rays.extend_from_slice(&[16.0 - 100.0 * dx, 16.0 - 100.0 * dy, 16.0, dx, dy, 0.0]);
            dirs.push(([16.0 - 100.0 * dx, 16.0 - 100.0 * dy, 16.0], [dx, dy, 0.0]));
        }
        // Oblique 3-D ray and a clean miss.
        let inv3 = 1.0 / (3.0f32).sqrt();
        rays.extend_from_slice(&[-20.0, -20.0, -20.0, inv3, inv3, inv3]);
        dirs.push(([-20.0, -20.0, -20.0], [inv3, inv3, inv3]));
        rays.extend_from_slice(&[-20.0, 500.0, 16.0, 1.0, 0.0, 0.0]);
        dirs.push(([-20.0, 500.0, 16.0], [1.0, 0.0, 0.0]));

        let step = 0.5f32;
        let mut gpu = vec![0.0f32; dirs.len()];
        projector
            .project_into(&device, &rays, step, &mut gpu)
            .expect("project");

        for (i, (o, d)) in dirs.iter().enumerate() {
            let cpu = Ray::try_from_direction(
                Point3::new(o[0], o[1], o[2]),
                Vector3::new(d[0], d[1], d[2]),
            )
            .and_then(|ray| forward_project_ray(&mu, &ray, step))
            .unwrap_or(0.0);
            // f32 GPU fma/mix vs CPU mul-add over ~10² steps: bound the relative
            // error at 1e-3 (≈ n·ε_f32 · condition margin), far below any physics
            // tolerance. Both sums are sequential per-ray (same order).
            assert!(
                (gpu[i] - cpu).abs() <= 1e-3 * (1.0 + cpu.abs()),
                "ray {i}: gpu {} vs cpu {cpu}",
                gpu[i]
            );
        }
        // The miss ray is exactly zero on both paths.
        assert_eq!(gpu[dirs.len() - 1], 0.0);
    }

    #[test]
    fn empty_batch_is_ok() {
        let Ok(device) = crate::default_device() else {
            return;
        };
        let projector = GpuProjector::new(&device, &mu_volume()).expect("upload");
        let mut out: [f32; 0] = [];
        projector
            .project_into(&device, &[], 0.5, &mut out)
            .expect("empty batch is a no-op");
    }
}
