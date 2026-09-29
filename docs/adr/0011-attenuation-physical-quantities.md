# ADR 0011: Type attenuation density inputs

- Status: Accepted

## Context

The CT attenuation-map solver accepted reference water density as a scalar
documented in g/cm³. The GPU attenuation mapper accepted both mass attenuation
coefficient and density as raw f32 values. Helios already uses Aequitas
`AreaPerMass` for mass attenuation, so these boundaries discarded physical
semantics before the CPU or GPU kernel.

## Decision

Accept Aequitas `MassDensity<T>` for the reference water density in
`helios_solver::attenuation_map`. Accept `AreaPerMass<f32>` and
`MassDensity<f32>` in `helios_gpu::GpuAttenuationMapper::new`. Convert to the
existing g/cm³ and cm²/g scalar representations exactly once at the CT
calibration or GPU-uniform boundary. Keep dense HU and linear-attenuation
volumes scalar because they are numerical storage contracts.

## Alternatives rejected

- Keep raw density/coefficient scalars: preserves a public unit gap.
- Add Helios-specific wrappers: duplicates Aequitas quantity ownership.
- Store quantities per voxel: changes dense-kernel representation without a
  field-descriptor requirement.

## Verification

- `helios_solver::attenuation_map` takes `MassDensity<T>`;
  `helios_gpu::GpuAttenuationMapper::new` takes `AreaPerMass<f32>` and
  `MassDensity<f32>`, and the existing attenuation-map analytical, f32,
  invalid-density, GPU differential, and end-to-end callers use the typed
  inputs.
- Rustfmt and source residue scans pass; the workspace compiles under the
  committed lock.
