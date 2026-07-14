//! Helios-domain example: `VoxelGrid` + dense `Volume<T>` construction,
//! trilinear sampling, and CT-number HU → μ-density preparation.
//!
//! Demonstrates that the universal Helios grid (`VoxelGrid::axis_aligned`)
//! owns the discrete-index ↔ world-coordinate map and that the dense
//! `Volume<T>` is generic over `T: Scalar` — built once, sampled either at
//! voxel index (fast, integer) or in world coordinates (trilinear, sub-voxel).
//!
//! Run with:  cargo run --example voxel_grid_construction -p helios-domain

use helios_domain::{Volume, VoxelGrid};
use helios_math::{Point3, Scalar};

/// Field affine in each axis: at voxel index `(i, j, k)` returns
/// `2·i + 3·j + 5·k + 7`, which is exactly reproducible by a trilinear
/// sample at any continuous world coordinate (affine fields reproduce exactly
/// under trilinear interpolation).
fn affine_value<T: Scalar>(idx: [usize; 3]) -> T {
    T::from_f64((2 * idx[0] + 3 * idx[1] + 5 * idx[2] + 7) as f64)
}

fn affine_volume<T: Scalar>() -> (VoxelGrid<T>, Volume<T>) {
    let grid = VoxelGrid::axis_aligned(
        [4, 5, 6],
        [T::from_f64(2.0), T::from_f64(3.0), T::from_f64(4.0)], // mm per axis
        Point3::new(T::from_f64(10.0), T::from_f64(20.0), T::from_f64(30.0)),
    )
    .expect("valid axis-aligned grid");
    let volume = Volume::from_shape_fn(grid, affine_value::<T>);
    (grid, volume)
}

/// Build an `f32` and an `f64` volume on identical topology and verify
/// trilinear sampling reproduces the affine field exactly under the analytic
/// oracle `2·ix + 3·iy + 5·iz + 7` (the same oracle the foundation tests use).
fn assert_trilinear_is_exact_affine<T: Scalar>() {
    let t = std::any::type_name::<T>();
    let (grid, vol) = affine_volume::<T>();
    assert_eq!(grid.dims(), [4, 5, 6]);
    assert_eq!(grid.num_voxels(), 120);

    // Voxel-center (sub-voxel midpoint, mid-grid) targets an interior coordinate.
    let checkpoint_xyz = Point3::new(T::from_f64(1.5), T::from_f64(2.25), T::from_f64(3.75));
    let actual = vol.sample_trilinear(checkpoint_xyz).unwrap();
    let expected = T::from_f64(2.0) * T::from_f64(1.5)
        + T::from_f64(3.0) * T::from_f64(2.25)
        + T::from_f64(5.0) * T::from_f64(3.75)
        + T::from_f64(7.0);
    let dx = (actual - expected).abs();
    assert!(
        dx <= T::from_f64(1e-9),
        "[{t}] |{actual:?} - {expected:?}| = {dx:?}"
    );

    // Integer vertex hits the voxel value exactly.
    let vert = Point3::new(T::from_f64(1.0), T::from_f64(2.0), T::from_f64(3.0));
    let actual = vol.sample_trilinear(vert).unwrap();
    let expected = affine_value::<T>([1, 2, 3]);
    assert_eq!(actual, expected);
    println!("[{t}] trilinear interpolation is exact on an affine field");
}

/// Out-of-grid and non-finite samples must surface `None` without
/// panicking — the boundary reject contract.
fn assert_out_of_bounds_is_none() {
    let (_, vol) = affine_volume::<f64>();
    assert_eq!(vol.sample_trilinear(Point3::new(-0.1, 0.0, 0.0)), None);
    assert_eq!(vol.sample_trilinear(Point3::new(3.5, 0.0, 0.0)), None); // nx-1 = 3
    assert_eq!(vol.sample_trilinear(Point3::new(f64::NAN, 0.0, 0.0)), None);
    assert_eq!(vol.get(99, 0, 0), None);
    println!("Out-of-grid and non-finite samples are rejected as None");
}

fn main() {
    assert_trilinear_is_exact_affine::<f64>();
    assert_trilinear_is_exact_affine::<f32>();
    assert_out_of_bounds_is_none();
    println!("VoxelGrid + Volume pair correctly under the trilinear sampling contract.");
}
