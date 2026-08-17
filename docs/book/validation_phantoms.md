# Chapter 23 — Reference Phantoms and Ground Truth

<!-- generated-figure-start -->
![Figure 35.1 — Reference Phantoms and Ground Truth](figures/ch35/fig01_35_reference_phantoms_and_ground_truth.svg)
*Figure 35.1 — Reference Phantoms and Ground Truth*
<!-- generated-figure-end -->

Helios includes synthetic test phantoms for validation and regression testing.

## Water Cylinder Phantom

The standard phantom for attenuation and dose validation:

```text
fn water_cylinder_phantom(nx: usize, spacing: f64) -> Volume<f64> {
    let c = (nx as f64 - 1.0) * spacing / 2.0;
    Volume::from_shape_fn(grid, move |idx| {
        let r = distance_from_axis(idx, c, spacing);
        if r <= 25.0 { 0.0 } else { -1000.0 }  // water cylinder in air
    })
}
```

## Bone Insert

A cortical bone insert (800 HU) within the water cylinder tests
heterogeneity correction in the collapsed-cone solver.

## Analytical Oracles

The helios-analysis crate provides:
- **Cylinder radon oracle**: exact Radon transform of a uniform cylinder
- **Exponential depth-dose**: pencil beam in homogeneous water

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Analytical Solutions and Regression Tests](validation_regression.md)
