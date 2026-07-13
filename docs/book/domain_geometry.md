# Voxel Grids and Volumetric Data

The core spatial data structure in Helios is the `VoxelGrid` + `Volume<T>` pair,
owned by `helios-domain`.

## VoxelGrid: The Geometry Contract

`VoxelGrid<T: Scalar>` encodes the discrete-index ↔ world-coordinate map for a
regular, axis-aligned 3-D grid:

```text
world[mm] = origin + (i·sx, j·sy, k·sz)
```

Construction validates all invariants upfront:

```rust
use helios_domain::VoxelGrid;
use helios_math::{Point3, Scalar};

let grid = VoxelGrid::axis_aligned(
    [256, 256, 64],             // nx × ny × nz voxels
    [f64::from_f64(1.0); 3],    // 1 mm isotropic spacing
    Point3::new(0.0, 0.0, 0.0), // origin at patient iso-centre
)?;

assert_eq!(grid.num_voxels(), 256 * 256 * 64);
```

## Volume\<T\>: Dense Scalar Storage

`Volume<T>` wraps a C-contiguous `leto::Array3<T>` backed by the same `VoxelGrid`.
All access is bounds-checked via `get(i, j, k)` → `Option<T>`, and hot kernels use
`as_slice()` for zero-copy iteration.

```rust
let volume = Volume::from_shape_fn(grid, |[i, j, k]| {
    T::from_f64((i + j + k) as f64)
});

// Trilinear sample at continuous voxel-index coords
let p = Point3::new(T::from_f64(1.5), T::from_f64(2.5), T::from_f64(0.5));
let val = volume.sample_trilinear(p).unwrap();
```

### Trilinear Sampling Contract

For an affine scalar field `f(i, j, k) = a·i + b·j + c·k + d`:
- Sampling at an integer vertex returns the **exact** stored value.
- Sampling at any continuous coordinate reproduces the affine field **exactly**
  (trilinear interpolation of an affine field is analytic).

## Generic over Scalar

Because `Volume<T: Scalar>` is generic over the eunomia `RealField`, the same
code runs at `f32` (GPU staging, memory-efficient) and `f64` (reference quality)
without separate types.

## Further Reading

- [Example: VoxelGrid and Volume Construction](examples/voxel_grid_construction.md)
- [`helios-domain` source](../../crates/helios-domain/src/)
