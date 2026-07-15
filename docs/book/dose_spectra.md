# Beam Hardening and Poly-Energetic Spectra

Real photon beams are poly-energetic: a 6 MV treatment beam contains photons
from ~0.1 MeV to 6 MeV. As the beam penetrates tissue, low-energy photons are
preferentially absorbed (beam hardening), shifting the mean energy upward.

## Spectral Components

Helios models a beam as a sum of SpectralComponent bins:

`ust
use helios_simulation::SpectralComponent;

let spectrum = vec![
    SpectralComponent { energy_mev: 1.0, weight: 0.35 },
    SpectralComponent { energy_mev: 3.0, weight: 0.45 },
    SpectralComponent { energy_mev: 6.0, weight: 0.20 },
];
`

## Beam-Hardened Transport

For each spectral component, the collapsed-cone kernel uses the component's
own attenuation coefficient μ(E), producing energy-dependent depth–dose
curves. The total dose is the energy-weighted sum:

`	ext
D(r) = Σ_k  w_k · D_k(r; μ_k)
`

## Clinical Significance

Ignoring beam hardening introduces depth-dose errors of 3–5% at 20 cm depth.
The poly-energetic model reduces this to < 1%.

## Further Reading

- [Collapsed-Cone Convolution](dose_collapsed_cone.md)
- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)