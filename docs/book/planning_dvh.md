# Chapter 15 — Dose-Volume Histograms

The DVH summarises the volumetric dose distribution of a target or organ:

`
ust
use helios_analysis::Dvh;

let dvh = Dvh::new(&dose, 200);  // 200 histogram bins
let d95 = dvh.d_percent(95.0);  // dose covering 95% of volume
let v20 = dvh.v_dose(20.0);     // volume receiving ≥ 20 Gy
`

## Cumulative DVH

The cumulative DVH V(d) gives the fraction of voxels receiving dose ≥ d:

`	ext
V(d) = |{r : D(r) ≥ d}| / |Volume|
`

## Clinical Metrics

| Metric | Meaning | Typical constraint |
|---|---|---|
| D95 | Dose to 95% of volume | ≥ 95% prescription dose (target) |
| D2 | Dose to hottest 2% | ≤ 107% prescription dose |
| V20 | Volume receiving 20 Gy | ≤ 35% (lung SBRT) |

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Gamma Index and Plan Verification](planning_gamma.md)