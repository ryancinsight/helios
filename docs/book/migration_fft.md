# Chapter 33 — Apollo: FFT and Spectral Methods

Helios uses spectral methods for collapsed-cone poly-energetic dose
calculation, regularized FBP reconstruction, and certain adjoint
reconstruction paths.  These migrate from `rustfft` (and `realfft` for
real-to-real) to **Apollo**'s forward FFT crate.  Apollo is *forward-only*
(no inverse FFT — every Atlas consumer can invert by conjugation) and is
wired for autodiff through [`coeus::Tensor`] so that a spectral solver
that has gradient information flows directly into optimization.

## The Apollo Surface

```rust
pub struct FftPlan {
    shape: Vec<usize>,
    inner: PlanInner,
}

pub trait FftEngine {
    fn plan(shape: &[usize]) -> FftPlan;
    fn forward_real(&self, x: &[f64], out: &mut [Complex64]);
    fn forward_complex(&self, x: &[Complex64], out: &mut [Complex64]);
    fn coeus_forward<T: CoeusScalar>(&self, x: &Tensor<T>) -> Tensor<T>;
}
```

A typical helios spectral solve becomes:

```rust
let plan = FftPlan::new(&[nx, ny, nz]);
let mut spectrum = vec![Complex64::zero(); nx * ny * nz];
plan.forward_real(&dose_field, &mut spectrum);

// Spectral kernel (poly-energetic mono-to-medium correction)
for (k, s) in spectrum.iter_mut().enumerate() {
    *s *= medium_response(k, energy_bins);
}

// IFFT via conjugation (Apollo is forward-only)
let recovered = plan.forward_complex(&spectrum, ...)
    .map(|c| c.conj())  // + normalization pass
    .collect();
```

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `rustfft::FftPlanner::new().plan_fft(N)` | `FftPlan::new(&[N])` |
| `realfft::RealFftPlanner::new()` | `plan.forward_real(x, out)` |
| manual real-to-complex + zero-pad | `forward_real` (one call) |
| `num_complex::Complex64` | `eunomia::ComplexField` impl |
| `rustfft::inverse` | forward + conjugate (Apollo) |

The legacy `Inverse` path disappears — Atlas handles inverse via forward
+ complex-conjugate, which is mathematically identical and unifies the
solver's autodiff path.

## Autodiff Bridge

When a solver needs a gradient (e.g. for adjoint optimization or
DVH-constrained treatment planning), Apollo returns a `Tensor<T>` from
[`coeus`]: the FFT preserves the autodiff tape so that a downstream loss
yields a backward through the FFT.  This is the critical feature that
motivates Apollo's design over `rustfft`.

## Helios-Specific FFT Use Cases

| Use case | Where |
|---|---|
| Poly-energetic collapsed-cone | `helios-solver/src/collapsed_cone.rs` |
| Regularized FBP ramp | `helios-imaging/src/fbp.rs` |
| MVCT forward + adjoint | `helios-imaging/src/mvct.rs` |
| Beam hardening correction | `helios-physics/src/spectra.rs` |

## Validation Examples

- [`collapsed_cone_3d`](examples/collapsed_cone_3d.md) — poly-energetic
  spectral dose engine.
- [`fbp_reconstruction`](examples/fbp_reconstruction.md) — FBP ramp
  filter.
- [`radon_sinogram`](examples/radon_sinogram.md) — Radon transform.

## Further Reading

- [`apollo` source](../../../apollo/)
- [Coeus: Tensors and Autodiff](migration_coeus.md)
- [`coeus` source](../../../coeus/)
