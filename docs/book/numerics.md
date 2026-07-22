# Chapter 3 — Scalar Fields and Numeric Abstractions

Helios is generic over the unomia::RealField numeric scalar, allowing
the same physics code to run at 32 (GPU/staging) or 64 (reference quality).

## The Scalar Hierarchy

`	ext
eunomia::NumericElement   ← integer types
         ↓
eunomia::FloatElement     ← f32, f64
         ↓
eunomia::RealField        ← algebraic field over ℝ
`

## Generic Physics

Helios domain objects use T: RealField:

`
ust
use helios_domain::Volume;
use helios_math::Scalar;

fn rms_dose<T: Scalar>(dose: &Volume<T>) -> T {
    let sum = dose.as_slice().iter().fold(T::zero(), |acc, &v| acc + v * v);
    (sum / T::from_usize(dose.num_voxels()).unwrap()).sqrt()
}
`

## Atlas Crate Integration

| Operation | Crate |
|---|---|
| Scalar traits | unomia |
| Array storage | leto::Array3<T> |
| SIMD dispatch | hermes-simd |
| GPU kernels | hephaestus-wgpu |

## Further Reading

- [Physics Domain Types and Safety Boundaries](foundations.md)
- [Memory and Allocation](memory.md)