# Chapter 27 — Eunomia: Numeric Trait Unification

Every Atlas crate depends on [`eunomia::RealField`] for scalar arithmetic and
on [`eunomia::ComplexField`] when complex spectral kernels appear.  Eunomia
is the **trait frontier** — a one-crate abstraction that the rest of the
stack builds on, replacing `num-traits`, `num-complex`, and the various
`num-derive` utilities.

## The RealField Trait

```rust
pub trait RealField:
    Copy
    + Send
    + Sync
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + PartialOrd
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_f64(x: f64) -> Self;
    fn to_f64(self) -> f64;
    fn sqrt(self) -> Self;
    fn abs(self) -> Self;
    fn min(self, other: Self) -> Self;
    fn max(self, other: Self) -> Self;
}
```

The trait is **deliberately small** — anything beyond the minimum lives in
extension traits.

## Companion Traits

```rust
pub trait FloatElement: RealField + ... { /* f32/f64-only rules */ }
pub trait IntElement:   ... { /* i32/i64-only rules */ }
pub trait ComplexField: ... { /* z = re + i*im rules */ }
```

These traits **carve the trait frontier in two**: `FloatElement` for
floating-point work (dose values, attenuation maps), `IntElement` for
indexing arithmetic, `ComplexField` for spectral methods.

## ZST and PhantomData Bounds

Atlas crates use ZSTs and `PhantomData` to encode capability requirements
**at the type level** with zero runtime cost:

```rust
pub struct DoseCalc<F: FloatElement, B: ComputeBackend> {
    inner: PhantomData<(F, B)>,
    grid:  VoxelGrid<F>,
}
```

`PhantomData<(F, B)>` is zero-sized — it compiles to nothing — but forces
the compiler to verify that `DoseCalc<f32, Cpu>` and `DoseCalc<f64, Wgpu>`
are **distinct types**.  This delivers **monomorphization without virtual
dispatch**.

## Migration From num-traits

| Legacy | Atlas |
|---|---|
| `T: num_traits::Float` | `T: eunomia::FloatElement` |
| `num_complex::Complex64` | `eunomia::ComplexField` impl |
| `from_f64` via `as` cast | `RealField::from_f64` |
| `as f32` / `as f64` at module boundaries | `RealField::from_f64` |

A typical helios port (e.g. for the dose engine):

```rust
// legacy: f64 only, with `as` casts
// atlas:
pub fn compute<F: FloatElement>(&self, grid: &VoxelGrid<F>) -> Result<DoseGrid<F>, HeliosError> {
    Ok(...)  // each F has a monomorphized kernel
}
```

## Validation Examples

- [`validate_foundation_units`](examples/validate_foundation_units.md) —
  validating newtypes consume `RealField::from_f64`.
- [Foundations chapter](foundations.md) — the `EnergyMeV`,
  `HounsfieldUnit`, `VoxelSpacingMm` types all parameterize over `RealField`.

## Further Reading

- [`eunomia` source](../../../eunomia/crates/)
- [Atlas Dependency Map](appendix_dependencies.md)
