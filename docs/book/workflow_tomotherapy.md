# Chapter 17 — TomoTherapy End-to-End Workflow

This chapter presents the complete helical TomoTherapy delivery and verification
pipeline — from CT imaging through dose calculation to plan QA.

## System Overview

```
CT acquisition → Attenuation map → MVCT reconstruction
                                           │
                      Helical MLC delivery ┤
                                           │
                                 Dose calculation (Collapsed-cone)
                                           │
                            DVH + Gamma-index verification
```

## CT and MVCT Imaging

The imaging path converts CT Hounsfield Units (HU) to linear attenuation
coefficients (μ) using mass-attenuation physics from `helios-physics`:

```rust
let mu = attenuation_map(&ct_volume, MassAttenuation::water());
let sinogram = parallel_beam_radon(&mu, n_angles);
let recon = filtered_back_projection(&sinogram, n_angles, nx);
```

## Helical Delivery Simulation

A helical MLC delivery is defined by a `LeafOpenTimeSinogram` (LOTS) encoding the
per-gantry-angle leaf opening times. `simulate_helical_delivery` produces the 3-D
terma (total energy released per unit mass) distribution:

```rust
let delivery = HelicalDelivery { sinogram: lots, beam_geometry: BeamGeometry::default() };
let terma = simulate_helical_delivery(&delivery, &mu_volume, &spectrum);
```

## Collapsed-Cone Dose

The `CollapsedCone` solver converts terma to dose using a poly-energetic
beam-hardening correction driven by the spectral components:

```rust
let dose = accumulate_delivered_dose_anisotropic(&terma, &mu, &CollapsedCone::default());
```

## Plan Quality Assurance

The `helios-analysis` crate provides:
- `Dvh` — Dose-Volume Histogram with D95/V20 lookup
- `gamma_index_3d` — per-voxel gamma metric (ΔD/δD)² + (Δr/δr)²
- `gamma_pass_rate` — fraction of voxels with γ < 1

```rust
let dvh = Dvh::new(&dose, 200);
let gamma = gamma_index_3d(&dose, &reference_dose, 0.03, 2.0, 5.0);
println!("Pass rate: {:.1}%", gamma_pass_rate(&gamma, 1.0) * 100.0);
```

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [`helios-simulation` source](../../crates/helios-simulation/src/)
- [`helios-analysis` source](../../crates/helios-analysis/src/)
