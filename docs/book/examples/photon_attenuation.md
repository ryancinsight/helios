# Example: Photon Attenuation Physics

**Crate**: `helios-physics`  
**Run**: `cargo run -p helios-physics --example photon_attenuation`  
**Source**: [`crates/helios-physics/examples/photon_attenuation.rs`](../../crates/helios-physics/examples/photon_attenuation.rs)

## What This Example Demonstrates

Exercises the core photon-interaction physics in `helios-physics`:

| Concept | API |
|---|---|
| Linear attenuation coefficient | `LinearAttenuation::new(mu)` |
| Beer–Lambert transmission | `mu.transmission(path_cm)` |
| Half-value layer | `mu.half_value_layer()` |
| Mass attenuation → linear | `MassAttenuation::to_linear(density)` |
| CT number → electron density | `relative_electron_density_from_hu(hu)` |
| CT number → mass density | `mass_density_from_hu(hu, water_rho)` |

## Key Code Snippet

```rust
use helios_physics::{LinearAttenuation, MassAttenuation, mass_density_from_hu};

// Water at 100 keV: μ ≈ 0.171 cm⁻¹
let mu = LinearAttenuation::new(0.171_f64).expect("valid μ");
let hvl = mu.half_value_layer().unwrap(); // ≈ 4.05 cm

// 10 cm path → transmitted fraction
let t = mu.transmission(10.0); // exp(-μ·x)

// CT-number calibration
let rho_rel = relative_electron_density_from_hu(0.0_f64); // water → 1.0
```

## Physics Background

The **Beer–Lambert law** describes narrow-beam attenuation:

```
I(x) / I₀ = exp(−μ · x)
```

where `μ` (cm⁻¹) is the linear attenuation coefficient. For a compound
material, `μ = (μ/ρ) · ρ` where `μ/ρ` (cm²/g) is the material-specific
mass attenuation coefficient (tabulated in NIST XCOM).

The **CT calibration** maps Hounsfield units to relative electron density via
`ρ_rel = max(0, 1 + HU/1000)` — the first-order linear model used as the
vendor-independent baseline before scanner-specific stoichiometric calibration.

## Book Chapter

[← Hounsfield Units and Attenuation Maps](../imaging_ct.md)
