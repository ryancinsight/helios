# ADR 0020: Advance the Eunomia and Moirai provider contracts

Status: Accepted  \
Date: 2026-11-15  \
Change class: [arch] [major]  \
Board item: [ATLAS-VERSION-COHERENCE-SWEEP-2026-09-30](../../backlog.md)

Revision 2026-11-16: the Leto upper bound moves from `<0.49.0` to
`<0.51.0` for the Leto 0.50.0 cohort. The lock holds Leto 0.50.0, which
CoeUs, RITK, Gaia's mesh crate, Hephaestus and Apollo now admit.

## Context

`helios-math` re-exports Eunomia's element vocabulary
(`RealField` as `Scalar`, `FloatElement`, `NumericElement`) and resolves
its parallel dispatch through Moirai's runtime and parallel executor.
The stack's provider cohort has moved past both lines:

- Eunomia 0.12.0 is the current provider head. Rust treats the
  `FloatElement` trait of Eunomia 0.11 and the same-named trait of
  Eunomia 0.12 as distinct nominal contracts, so a consumer that names
  the re-exported trait while depending on Eunomia 0.11 directly
  resolves two distinct traits for the same name, and the mixture stops
  compiling.
- Moirai 0.12.0 removes the required `StateWord::get_mut` method.
  Helios implements no `StateWord` and calls no removed method, so no
  Helios source changes with the move.

The requirement also gates the rest of the stack. The providers' current
default heads require Moirai 0.12 (`apollo-fft` 0.29.5 and
`coeus-core` 0.12.0 name `moirai-runtime ^0.12.0`), so every consumer
that resolves those providers with a `^0.11.0` requirement cannot
re-resolve onto the current provider heads at all.

Two alternatives were considered:

1. Stay on Eunomia 0.11 and Moirai 0.11 and keep the current lock. The
   provider cohort then splits across two Eunomia identities and two
   Moirai identities inside one workspace, and the consumers that
   resolve Helios inherit the same split.
2. Require `>=0.10.0, <0.13.0` for Eunomia and let the lock choose. The
   range admits both provider lines, but the lock selects the newest
   compatible revision anyway, so the re-exported traits move exactly
   as they do under the selected option. The wider range only hides the
   move from the requirement.

## Decision

Move the Eunomia requirement to 0.12.0, the Moirai runtime and
Moirai-parallel requirements to 0.12.0, and widen the Leto requirement
to `>=0.46.4, <0.51.0` so the 0.48, 0.49 and 0.50 lines resolve. Regenerate the lock
outside the Atlas dependency overlay, so every provider row names one
revision per provider repository and the `source = "git+..."` rows
survive.

No Helios source changes with the move. The removed provider method is
unused in this repository.

## Consequences and migration

- Code generic over the Eunomia 0.11 `RealField` trait no longer
  satisfies Helios's public bounds. A consumer updates its Eunomia bound
  and imports to 0.12 and selects the same provider cohort before it
  uses Helios's public generic types. Consumers already using concrete
  scalar types keep the same Helios API after their dependency graph
  selects the new cohort.
- The public contract change is a [major] change class; the release
  tooling assigns the version at release.
- No compatibility trait, forwarding wrapper, or duplicate element
  vocabulary is added: two Eunomia identities under one name would keep
  the split the provider owns.
