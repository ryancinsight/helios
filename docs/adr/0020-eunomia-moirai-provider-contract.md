# ADR 0020: Advance the Eunomia and Moirai provider contracts

Status: Accepted  \
Date: 2026-11-15  \
Change class: [arch] [major]  \
Delivered by: PR #199

Revision 2026-11-23: the Decision no longer states a Leto requirement.
PR #214 replaced it with `0.51.0` and PR #219 moved it to `0.52.0`, so
the range this record gave was stale. [ADR 0022] records the Leto 0.52
decision and owns the requirement. The Decision also names the Eunomia
0.12 line (`Cargo.toml` requires 0.12.1) and the registry lock source
that [ADR 0022] records, in place of the `git+` rows. Context dates its
provider-head claims. Consequences now name the Leto types Helios
re-exports. PR #227 carries this revision.

## Context

`helios-math` re-exports Eunomia's element vocabulary
(`RealField` as `Scalar`, `FloatElement`, `NumericElement`) and resolves
its parallel dispatch through Moirai's runtime and parallel executor.
At the decision date the stack's provider cohort had moved past both
lines:

- Eunomia 0.12.0 was the current provider head. Rust treats the
  `FloatElement` trait of Eunomia 0.11 and the same-named trait of
  Eunomia 0.12 as distinct nominal contracts, so a consumer that names
  the re-exported trait while depending on Eunomia 0.11 directly
  resolves two distinct traits for the same name, and the mixture stops
  compiling.
- Moirai 0.12.0 removes the required `StateWord::get_mut` method.
  Helios implements no `StateWord` and calls no removed method, so no
  Helios source changes with the move.

The requirement also gates the rest of the stack. The providers' default
heads then required Moirai 0.12 (`apollo-fft` 0.29.5 and `coeus-core`
0.12.0 named `moirai-runtime ^0.12.0`), so every consumer that resolved
those providers with a `^0.11.0` requirement could not re-resolve onto
the provider heads at all.

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

Move the Eunomia requirement to the 0.12 line (`Cargo.toml` requires
`0.12.1`) and the Moirai runtime and Moirai-parallel requirements to
0.12.0. [ADR 0022] owns the Leto requirement. Regenerate the lock
outside the Atlas dependency overlay, so it resolves every provider from
the published registry source that [ADR 0022] records, with one version
per provider crate.

No Helios source changes with the move. The removed provider method is
unused in this repository.

## Consequences and migration

- Code generic over the Eunomia 0.11 `RealField` trait no longer
  satisfies Helios's public bounds. A consumer updates its Eunomia bound
  and imports to 0.12 and selects the same provider cohort before it
  uses Helios's public generic types. Consumers already using concrete
  scalar types keep the same Helios API after their dependency graph
  selects the new cohort.
- Helios re-exports Leto's geometry vocabulary and returns it from a
  public method, so the Leto line is part of the public contract beside
  the Eunomia line. `crates/helios-math/src/lib.rs:66-69` re-exports
  `Isometry3`, `Point2`, `Point3`, `RotationBasisError`, `Translation3`,
  `UnitQuaternion`, `UnitVector3` and `Vector3` from `leto::geometry`
  under the `geometry` feature, and
  `crates/helios-domain/src/grid.rs:162` returns `Isometry3<T>` from
  `VoxelGrid::pose`. Each Leto line gives these types distinct Rust
  identities, so a consumer selects the Leto line that [ADR 0022]
  records.
- The public contract change is a [major] change class; the release
  tooling assigns the version at release.
- No compatibility trait, forwarding wrapper, or duplicate element
  vocabulary is added: two Eunomia identities under one name would keep
  the split the provider owns.

[ADR 0022]: 0022-published-physics-provider-boundary.md
