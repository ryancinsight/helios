# ADR 0010: Type Compton photon energy

- Status: Accepted
- Date: 2026-08-06
- Scope: `helios-physics` Compton kernels and the `helios-python` adapter

## Context

The public Klein–Nishina and Compton mass-coefficient functions accepted a
scalar documented as MeV. The scalar carried no energy dimension, so a caller
could pass joules, electronvolts, or an arbitrary numeric value without a type
error. Helios already uses Aequitas quantities for dose, fluence, and geometry.

## Decision

Accept `aequitas::systems::si::quantities::Energy<T>` at every public Compton
photon-energy boundary. Convert to `MegaElectronVolt` exactly once at the
dimensionless Klein–Nishina kernel boundary. Keep electron density, angular
cosines, cross-sections, and transfer fractions in their existing scalar or
typed result contracts.

The Python binding retains its explicit MeV scalar input because Python has no
compile-time quantity type; it validates the scalar and constructs an Aequitas
`Energy<f64>` before calling the Rust core. Examples and Rust tests construct
energy values through the Aequitas unit API.

## Alternatives rejected

- Keep a raw MeV scalar: preserves the dimensional gap.
- Add a Helios `EnergyMeV` wrapper around the public physics API: duplicates
  Aequitas ownership and prevents callers from selecting another energy unit.
- Convert every scalar at each formula: spreads unit-boundary logic through the
  numerical kernel instead of keeping one explicit conversion seam.

## Verification

- The Compton unit-equivalence test compares 1 MeV and 1,000,000 eV inputs.
- Existing analytical, monotonicity, f32, property, and NIST coefficient tests
  remain the behavioral oracles for the unchanged formulas.
- The Python and runnable example call sites construct the typed core input.
- The current Atlas metadata graph resolves `coeus-autograd` from
  `repos/coeus/crates/coeus-autograd`; no member-local patch or compatibility
  path is retained. Full warning-denied Clippy, 286/286 workspace Nextest,
  workspace doctest, warning-denied Rustdoc, and formatting gates pass.
