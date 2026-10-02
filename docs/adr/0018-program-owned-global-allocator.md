# ADR 0018: Programs own the global allocator; libraries depend only on what they use

- Status: Accepted
- Date: 2026-09-26
- Scope: `helios-allocator`, every example, benchmark, and binary; the
  Mnemosyne dependency edges of `helios-quantities`, `helios-planning`,
  `helios-python`, and `xtask`
- Delivered by: PR #112

## Context

PR #112 adopted Mnemosyne as the Helios allocator by making
`mnemosyne-memory` a default-on dependency of `helios-quantities` behind a feature
named `mnemosyne-memory`, with a `helios_quantities::install_global_allocator!`
macro. Review found four consequences:

- every library consumer of `helios-quantities` compiled the allocator stack:
  `cargo tree -e normal -p helios-quantities` grew from 24 to 46 packages
  (72 lines), and `helios-python` linked it without installing it;
- the feature was named after a stack member and existed only to switch a
  dependency off, while the stack's standards admit a capability into the one
  build or delete it;
- `helios-planning` reached `AlignedVec` through `mnemosyne-arena` while
  `helios_quantities::memory` re-exported the same type, so one type had two paths
  and the "single naming site" claim was false;
- `xtask` installed the allocator although `bench-replicated` measures a
  separately spawned `cargo bench` process, and `ramp_filter` did not install
  it at all.

## Decision

1. The global allocator is a program concern. A dedicated crate,
   `helios-allocator`, depends on `mnemosyne-memory` and exports one macro,
   `install_global_allocator!()`, expanding to a `#[global_allocator]` static.
   Every example and criterion bench invokes it once at its crate root; the
   owning library crate lists `helios-allocator` under `[dev-dependencies]`, so
   the edge never reaches that library's consumers.
2. Library crates depend only on the capability they use. `helios-quantities`
   carries no Mnemosyne dependency and no memory module. `helios-planning`
   depends on `mnemosyne-arena` directly, because its dose-influence matrix
   stores entries in `AlignedVec`; routing that type through `helios-quantities`
   would put `mnemosyne-arena` under every `helios-quantities` consumer to serve one
   crate.
3. No feature gates any of this: each crate's graph holds exactly the
   dependencies its code uses.
4. `helios-python` installs no allocator — it is a `cdylib` inside CPython —
   and, with (2), links no allocator crate; it keeps `mnemosyne-arena` only
   through `helios-planning`'s storage.
5. `xtask` installs no allocator: nothing it measures runs in its process.

## Rejected alternatives

- **Default feature on `helios-quantities`** (PR #112 as opened): the cost lands on
  every library consumer, and disabling it needs a member-named feature.
- **Re-export `AlignedVec` from `helios-quantities`**: one path per type, but it
  adds `mnemosyne-arena` to all 24 packages' consumers for one user.
- **A two-line `#[global_allocator]` static in each program** naming
  `mnemosyne` directly: 25 copies of the choice instead of one.

## Consequences

- Changing the allocator edits `helios-allocator` alone.
- A new example or bench must invoke the macro; the doc chapter and
  `ARCHITECTURE.md` state that rule.
- Mnemosyne's `GlobalAlloc` contract (alignment, zero-size, `realloc`) is
  tested in Mnemosyne (`crates/mnemosyne/tests/global_alloc_tests/`); the
  `helios-allocator` doctest runs a program under the installed allocator.

## Overturning evidence

A library that must itself allocate through Mnemosyne-specific APIs (arena
placement, scratch pools) takes a direct dependency on the Mnemosyne crate
providing that API, under (2) — it does not reopen (1).
