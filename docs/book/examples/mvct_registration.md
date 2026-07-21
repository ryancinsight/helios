# Example: IGRT Setup Correction via Translation Registration

**Crate**: `helios-imaging`  
**Run**: `cargo run -p helios-imaging --example mvct_registration`  
**Source**: [`crates/helios-imaging/examples/mvct_registration.rs`](../../crates/helios-imaging/examples/mvct_registration.rs)

## What This Example Demonstrates

Integer-voxel translation registration for image-guided radiation therapy (IGRT)
on a 32×32×32 synthetic water/bone phantom with a known applied setup error.

| Stage | API | Outcome |
|---|---|---|
| Planning CT | `Volume::from_shape_fn` | 10×10×10 bone insert at centre |
| Shifted daily MVCT | Manual index offset | +3, −2, +1 voxel setup error |
| Registration | `register_translation(&reference, &daily, max_shift)` | Detected shift exact |
| Couch correction | Apply `−detected_shift` | RMSE < 1e-6 cm⁻¹ |

## Key Code Snippet

```rust
use helios_imaging::register_translation;

// Max search radius ±5 voxels per axis (covers ≤15 mm at 3 mm resolution)
let detected_shift = register_translation(&reference, &daily, [5, 5, 5]);
// → [3, -2, 1] exactly matching the applied error

// Convert to mm for couch shift table
let couch_mm: Vec<f64> = detected_shift.iter().map(|&s| s as f64 * voxel_mm).collect();
```

## Physics Background

Before radiotherapy delivery the patient position is verified with a daily MVCT
and aligned to the planning CT. The whole-voxel translation registrar minimizes
the mean squared intensity difference over the overlap region:

```
cost(s) = mean_v ( daily(v) − reference(v − s) )²
```

For a uniform phantom translated by integer voxels the minimum is exactly zero
at the true shift — the deterministic, analytically verifiable base case.
Sub-voxel shifts require interpolation (helios H-044b); rotation and deformable
registration scale up through `ritk` mutual-information.

## Book Chapter

[← MVCT and Correction Workflows](../imaging_mvct.md)
