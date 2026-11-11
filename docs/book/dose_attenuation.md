# Chapter 9 — Mass Attenuation and Photon Cross Sections

<!-- generated-figure-start -->
![Figure 9.1 — Mass Attenuation and Photon Cross Sections](figures/ch09/fig01_9_mass_attenuation_and_photon_cross_sections.svg)
*Figure 9.1 — Mass Attenuation and Photon Cross Sections*
<!-- generated-figure-end -->

Hyperion provides the MassAttenuation coefficient models, and
helios-solver maps Hounsfield-unit volumes to linear attenuation
coefficients (μ, cm⁻¹):

```text
use aequitas::systems::si::{
    quantities::{AreaPerMass, MassDensity as DensityQuantity},
    units::{GramPerCubicCentimeter, SquareCentimeterPerGram},
};
use helios_solver::attenuation_map;
use hyperion::coefficient::MassAttenuation;

let mu_over_rho = MassAttenuation::new(AreaPerMass::from_unit::<SquareCentimeterPerGram>(0.0636))?;
let water = DensityQuantity::from_unit::<GramPerCubicCentimeter>(1.0);
let mu = attenuation_map(&ct_volume, mu_over_rho, water)?;
```

## Water Calibration

For photon beams in the MV range, the linear attenuation coefficient is
derived from the HU value and the beam energy:

```text
μ(HU) = μ_water · (1 + HU/1000)   for soft tissue
μ(HU) = μ_cortical_bone            for HU > 700
```

Typical values at 6 MV:
| Material | HU | μ (cm⁻¹) |
|---|---|---|
| Air | −1000 | ~0 |
| Water | 0 | 0.060 |
| Cortical bone | 800 | 0.108 |

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Terma and Energy Deposition](dose_terma.md)
