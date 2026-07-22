# Example: VoxelGrid and Volume Construction

**Crate**: `helios-domain`  
**Run**: `cargo run -p helios-domain --example voxel_grid_construction`  
**Source**: [`crates/helios-domain/examples/voxel_grid_construction.rs`](../../../crates/helios-domain/examples/voxel_grid_construction.rs)

## What This Example Demonstrates

This example builds a 4×5×6 voxel grid with non-uniform spacing (2 mm, 3 mm, 4 mm
per axis) and fills it with the affine scalar field `f(i,j,k) = 2i + 3j + 5k + 7`.
It then verifies the trilinear sampling contract:

1. **Affine exactness** — trilinear interpolation of an affine field reproduces the
   analytic value exactly (to 1e-9 absolute error for `f64`, exact for integer indices).
2. **Out-of-bounds rejection** — coordinates outside `[0, nx−1]×[0, ny−1]×[0, nz−1]`
   and non-finite coordinates return `None` without panicking.
3. **Generic over `Scalar`** — the same test runs with `f64` and `f32`.

## Key Code Snippet

```rust
use helios_domain::{Volume, VoxelGrid};
use helios_math::{Point3, Scalar};

fn affine_volume<T: Scalar>() -> (VoxelGrid<T>, Volume<T>) {
    let grid = VoxelGrid::axis_aligned(
        [4, 5, 6],
        [T::from_f64(2.0), T::from_f64(3.0), T::from_f64(4.0)],
        Point3::new(T::from_f64(10.0), T::from_f64(20.0), T::from_f64(30.0)),
    ).expect("valid axis-aligned grid");
    let volume = Volume::from_shape_fn(grid, |[i, j, k]| {
        T::from_f64((2*i + 3*j + 5*k + 7) as f64)
    });
    (grid, volume)
}
```

## Trilinear Contract Proof

For an affine field `f(x) = a·x + b`, trilinear interpolation at continuous coordinate
`x_c ∈ [x₀, x₁]` gives:
```
lerp(f(x₀), f(x₁), t) = (1-t)·(a·x₀+b) + t·(a·x₁+b) = a·((1-t)x₀ + t·x₁) + b = f(x_c)
```

The same argument applies in 3D (each lerp axis preserves affine exactness), so
`sample_trilinear` on an affine volume is analytically exact — any error comes only
from floating-point rounding.

## Book Chapter

[← Voxel Grids and Volumetric Data](../domain_geometry.md)
