# ADR 0002: Proteus material-density boundary

- Status: Accepted
- Date: 2026-07-20
- Class: [arch] [major]

## Context

`MassAttenuation::to_linear` duplicated finite, non-negative mass-density
validation even though density is a material property shared with acoustic and
fluid solvers. Helios owns photon-interaction and CT-calibration laws, not the
cross-domain material-property validity boundary.

## Decision

- Accept Proteus `MassDensity<T>` at the mass-to-linear attenuation boundary.
- Convert its Aequitas quantity to g/cm3 exactly once for the Helios-owned
  `mu = (mu/rho) rho` law.
- Keep CT-number calibration and mass-attenuation coefficients in Helios.
- Update callers in the same change; retain no raw-density overload.

## Alternatives rejected

- Keep raw generic scalars and repeat validation: rejected because the same
  material invariant would retain multiple owners.
- Move attenuation into Proteus: rejected because photon interaction is not a
  shared constitutive material law in the current consumer set.
- Add a forwarding overload: rejected because it would preserve the obsolete
  unvalidated contract.

## Consequences

The public method is intentionally breaking. Invalid density cannot reach the
attenuation law through safe construction. Unit conversion remains typed
through Aequitas, and the transparent Proteus newtype adds no storage or
dispatch cost.
