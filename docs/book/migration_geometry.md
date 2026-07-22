# Chapter 29 — Leto: Geometry — VoxelGrid, MLC, Beam Isometries

Helios geometry — voxel grids, multi-leaf collimator (MLC) leaf positions,
beam isometries, helical delivery carousels — migrates from
`nalgebra::Point3`/`Vector3`/`Isometry3` to the matching **Leto** types.
Leto's geometry types share the same in-memory layout and trait surface as
their `nalgebra` analogues, but encode compile-time distinctions between
scalar types (`f32` vs `f64`), precision wrappers, and backend (CPU/GPU)
via `PhantomData` so that callers cannot accidentally cross wires.

## The Core Types

```rust
pub struct Point3<T: RealField>       { data: [T; 3], _marker: PhantomData<()> }
pub struct Vector3<T: RealField>      { data: [T; 3], _marker: PhantomData<()> }
pub struct Rotation3<T: RealField>    { data: Quaternion<T>, _marker: PhantomData<()> }
pub struct Translation3<T: RealField> { data: Vector3<T>, _marker: PhantomData<()> }
pub struct Isometry3<T: RealField>    { data: Rotation3<T>, translation: Vector3<T> }
```

Each type:

- Stores three scalar fields contiguously (`size_of::<Point3<f64>>() == 24`).
- Carries a ZST `PhantomData` to force **monomorphization** per `(F,
  backend)` pair.
- Implements `Copy` when `T: Copy` (true for `f32`/`f64`), so passing
  geometry into a kernel is **zero-clone**.

## Migration From nalgebra

| Legacy | Atlas |
|---|---|
| `nalgebra::Point3<f64>` | `leto::Point3<f64>` |
| `nalgebra::Vector3<f64>` | `leto::Vector3<f64>` |
| `nalgebra::Isometry3<f64>` | `leto::Isometry3<f64>` |
| `Matrix3::face_normal(...)` | `Vector3::cross(a, b).normalize()` |

## Helios-Specific Geometry: VoxelGrid

```rust
pub struct VoxelGrid<F: FloatElement> {
    shape: [u32; 3],            // (nx, ny, nz)
    spacing_mm: [F; 3],         // voxel pitch [mm]
    origin_mm: Point3<F>,       // patient-coordinate origin
}

impl<F: FloatElement> VoxelGrid<F> {
    pub fn axis_aligned(shape: [u32; 3], spacing_mm: [F; 3], origin_mm: Point3<F>) -> Self;
    pub fn index_to_point(&self, i: u32, j: u32, k: u32) -> Point3<F>;
}
```

VoxelGrid replaces ad-hoc `Vec<f64>` × `(nx, ny, nz)` triples.  The
typed `origin_mm` and `spacing_mm` cross-check units at construction,
preventing the silent mm-vs-cm bugs that dominated the legacy code base.

## Helios-Specific Geometry: MLC

```rust
pub struct Mlc<F: FloatElement> {
    leaves: Vec<MlcLeaf<F>>,    // 1 leaf per bank per row
    bank_count: u8,
    rows: u16,
}

pub struct MlcLeaf<F: FloatElement> {
    position_mm: F,             // signed: +ve extension, -ve retraction
    width_mm:    F,
}
```

MLC replaces a stack of `Vec<(f64, f64)>` extent tuples plus a runtime
"is this a Varian-style or Elekta-style" enum.  The typed `bank_count` and
`rows` prevent off-by-one leaf sequencing.

## Beam Isometries

```rust
pub struct BeamIsometry<F: FloatElement> {
    gantry_to_patient: Isometry3<F>,    // gantry rotation
    collimator_to_patient: Isometry3<F>,// collimator rotation
    couch_to_patient: Isometry3<F>,     // couch offset
}
```

Beam-isometry chains save the legacy code's per-component matrix products
(3 × 3 = 9 multiplications per call site).  One `Isometry3<T>::compose`
chain returns the right transformation with type-monomorphized math.

## Validation Examples

- [`voxel_grid_construction`](examples/voxel_grid_construction.md) —
  VoxelGrid construction and index-to-point lookup.
- [Domain Geometry chapter](domain_geometry.md) — VoxelGrid, Volume,
  storage abstraction.
- [Treatment Planning chapters](planning_mlc.md), (planning_helical.md).

## Further Reading

- [`helios-domain` source](../../crates/helios-domain/src/)
- [`leto` geometry module](../../../leto/crates/)
- [Leto: Arrays](migration_arrays.md)
