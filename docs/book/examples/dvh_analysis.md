# Example: DVH Analysis

**Crate**: `helios-analysis`  
**Run**: `cargo run -p helios-analysis --example dvh_analysis`  
**Source**: [`crates/helios-analysis/examples/dvh_analysis.rs`](../../crates/helios-analysis/examples/dvh_analysis.rs)

## What This Example Demonstrates

Builds a synthetic Gaussian dose distribution and computes clinical DVH
evaluation metrics used in treatment plan review:

| Metric | API | Clinical Meaning |
|---|---|---|
| D₉₅ | `dvh.dose_at_volume_fraction(0.95)` | 95% of volume receives ≥ D₉₅ |
| V₉₅% | `dvh.volume_fraction_at_dose(Rx·0.95)` | fraction receiving ≥ 95% Rx |
| Homogeneity Index | `dvh.homogeneity_index()` | ICRU-83 (D₂−D₉₈)/D₅₀ |
| Masked DVH | `Dvh::from_volume_masked(dose, mask)` | structure-specific DVH |

## Key Code Snippet

```rust
use helios_analysis::Dvh;

let dvh = Dvh::from_volume(&dose_volume);
println!("D₉₅ = {:.2} Gy", dvh.dose_at_volume_fraction(0.95));
println!("HI  = {:.4}",    dvh.homogeneity_index());

// PTV-masked DVH
let ptv_dvh = Dvh::from_volume_masked(&dose, |[i, j, k]| ptv_mask[i][j][k]);
```

## Clinical Context

The **dose-volume histogram** is the primary clinical plan evaluation tool.
Key metrics per ICRU-83 and AAPM TG-119:

- **Coverage**: `D₉₅ ≥ 95% Rx` (PTV coverage gate)
- **Hotspot**: `D₂ ≤ 107% Rx` (maximum dose constraint)
- **Homogeneity Index**: `HI = (D₂ − D₉₈) / D₅₀ ≤ 0.10` (uniform dose ideal)
- **Structure-specific** DVHs are computed with `from_volume_masked` for each
  PTV and OAR using the RT-struct rasterization mask.

## Book Chapter

[← Dose-Volume Histograms](../planning_dvh.md)
