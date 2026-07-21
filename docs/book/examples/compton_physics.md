# Example: Compton Scattering Physics

> **Source:** `crates/helios-physics/examples/compton_physics.rs`
>
> **Run:** `cargo run -p helios-physics --example compton_physics`

## Overview

Demonstrates the Klein–Nishina cross-section model that underpins MV photon transport in `helios-physics`. At megavoltage beam energies Compton scattering dominates (photoelectric and pair-production are secondary), so the Klein–Nishina formula is the central physics primitive.

| Verification | Expected result |
|---|---|
| KN → Thomson at 1 keV | ratio ≈ 1.0 (within 1 %) |
| Cross section monotonically decreasing | σ_KN(E₁) > σ_KN(E₂) for E₁ < E₂ |
| Energy-transfer fraction > 0.5 at 6 MeV | `f_tr(6 MeV) ≈ 0.644` |

## Key APIs

```rust
use helios_physics::{
    thomson_cross_section, klein_nishina_cross_section,
    compton_mean_energy_transfer_fraction, compton_mass_attenuation,
    compton_energy_transfer_cross_section, electrons_per_gram,
};

let sigma_t: f64 = thomson_cross_section();             // 6.65e-29 m²/e⁻
let sigma_kn = klein_nishina_cross_section(6.0_f64);    // 6 MeV
let f_tr = compton_mean_energy_transfer_fraction(6.0);  // ≈ 0.64
let e_per_g = electrons_per_gram(0.5551_f64);           // water Z/A
let mu_rho = compton_mass_attenuation(0.1_f64, e_per_g); // cm²/g at 100 keV
```

## Physics Background

The Klein–Nishina total cross-section per electron is:
```text
σ_KN = 2πr_e² { (1+α)/α² [ 2(1+α)/(1+2α) - ln(1+2α)/α ]
                + ln(1+2α)/(2α) - (1+3α)/(1+2α)² }
```
where `α = E / m_e c²` (photon energy in electron-rest-mass units, `m_e c² = 0.511 MeV`) and `r_e = 2.818 × 10⁻¹⁵ m` is the classical electron radius.

As `α → 0` this reduces to the Thomson limit `σ_T = (8/3)π r_e² ≈ 6.652 × 10⁻²⁹ m²`.

## Book Chapter

[← Mass Attenuation and Photon Cross Sections](../dose_attenuation.md)
