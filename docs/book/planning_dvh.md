# Chapter 15 — Dose-Volume Histograms

The DVH summarises the volumetric dose distribution of a target or organ:

```rust
use aequitas::systems::si::{quantities::AbsorbedDose, units::Gray};
use helios_analysis::Dvh;

let dvh = Dvh::from_volume(&dose);
let d95 = dvh.dose_at_volume_fraction(0.95); // typed dose
let v20 = dvh.volume_fraction_at_dose(AbsorbedDose::from_unit::<Gray>(20.0));
```

## Cumulative DVH

The cumulative DVH V(d) gives the fraction of voxels receiving dose ≥ d:

```text
V(d) = |{r : D(r) ≥ d}| / |Volume|
```

## Clinical Metrics

| Metric | Meaning | Typical constraint |
|---|---|---|
| D95 | Dose to 95% of volume | ≥ 95% prescription dose (target) |
| D2 | Dose to hottest 2% | ≤ 107% prescription dose |
| V20 | Volume receiving 20 Gy | ≤ 35% (lung SBRT) |

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Gamma Index and Plan Verification](planning_gamma.md)
