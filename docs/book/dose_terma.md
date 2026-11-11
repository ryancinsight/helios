# Chapter 10 — Terma and Energy Deposition

<!-- generated-figure-start -->
![Figure 10.1 — Terma and Energy Deposition](figures/ch10/fig01_10_terma_and_energy_deposition.svg)
*Figure 10.1 — Terma and Energy Deposition*
<!-- generated-figure-end -->

**Terma** (Total Energy Released per unit MAss) is the primary quantity
from which dose is convolved:

```text
T(r) = Φ(r) · (μ_en/ρ)(r)
```

where Φ is the fluence and μ_en/ρ is the energy-absorption coefficient.

```text
use aequitas::systems::si::{
    quantities::{Angle, Dimensionless, Length, Time},
    units::{Millimeter, Radian, Second},
};
use helios_domain::{HelicalDelivery, LeafOpenTimeSinogram, MlcModel};
use helios_simulation::simulate_helical_delivery;

let delivery = HelicalDelivery::new(
    51,
    Length::from_unit::<Millimeter>(20.0),
    Dimensionless::from_base(0.3),
    Time::from_unit::<Second>(10.0),
    Angle::from_unit::<Radian>(0.0),
    Length::from_unit::<Millimeter>(0.0),
)?;
let lot = LeafOpenTimeSinogram::from_fractions(projections, leaves, fractions)?;
let mlc = MlcModel::new(0.01, 0.05)?;
let frames = simulate_helical_delivery(&delivery, &lot, &mlc);
```

## Helical Delivery Terma

For a TomoTherapy beam rotating helically, the terma accumulates from
all gantry angles and couch positions. The implementation in
helios-simulation sums fan-beam contributions at each delivery step.

## Physical Units

- Terma: Gy = J/kg
- Fluence Φ: photons/m²
- μ_en/ρ: m²/kg

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Collapsed-Cone Convolution](dose_collapsed_cone.md)
