# Chapter 11 — Collapsed-Cone Convolution

<!-- generated-figure-start -->
![Figure 11.1 — Collapsed-Cone Convolution](figures/ch11/fig01_11_collapsed_cone_convolution.svg)
*Figure 11.1 — Collapsed-Cone Convolution*
<!-- generated-figure-end -->

The CollapsedCone solver in helios-simulation converts terma to dose
by convolving with a poly-energetic pencil-beam kernel:

```text
use aequitas::systems::si::{quantities::Length, units::Millimeter};
use helios_simulation::accumulate_delivered_dose_anisotropic;

let dose = accumulate_delivered_dose_anisotropic(
    &frames,
    &mu,
    geometry,
    Length::from_unit::<Millimeter>(3.0),
    Length::from_unit::<Millimeter>(0.25),
    &cone,
)?;
```

## Algorithm

1. **Collapse** the pencil-beam point-spread function onto N_cone discrete
   angular cones (default: 48 cones)
2. **Ray-march** each cone direction through the attenuation map
3. **Accumulate** dose from exponential transport along each ray

```text
D(r) = Σ_k  T(r) · A_k · exp(−μ_bar · d_k(r))
```

## Beam Hardening

For poly-energetic beams, each spectral component is transported
with its own range kernel and energy-fluence weight, then summed:

```text
use aequitas::systems::si::{
    quantities::{Dimensionless, Length},
    units::Centimeter,
};
use helios_simulation::SpectralComponent;

let spectrum = vec![SpectralComponent {
    range_up: Length::from_unit::<Centimeter>(0.05),
    range_down: Length::from_unit::<Centimeter>(0.90),
    weight: Dimensionless::from_base(1.0),
}];
```

## Performance

The anisotropic variant uses moirai for parallel cone-direction dispatch.

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Beam Hardening and Poly-Energetic Spectra](dose_spectra.md)
