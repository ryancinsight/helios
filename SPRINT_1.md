# Sprint 1 — Foundation

**Goal:** Establish the Helios integrator workspace on the Atlas stack — skeleton,
foundation crate, and the architectural model — following the kwavers/cfdrs process.
Then begin `helios-core` → `helios-math` → `helios-domain` (ritk/gaia) build-out.

**Phase:** 1 (Foundation, 0–10%): 100% audit, planning, gap analysis. Transitioning
into Phase 2 (Execution) as the first crate lands.

## Decisions (ADR-lite)

1. **Repository role.** Helios is an *integrator* repo in the Atlas multi-repo stack
   (sibling to kwavers), a Cargo workspace consuming Atlas crates as remote git
   dependencies. Each Atlas repo carries its own `[workspace]`; no path scaffolding.
   - *Alternative rejected:* vendoring Atlas crates as path deps — breaks independent
     versioning and the co-evolution protocol.
2. **Edition 2021 / resolver 2.** Explicit goal directive; matches kwavers. Overrides
   the standards' edition-2024 default (recorded).
3. **Layering.** Strict unidirectional graph, `helios-core` innermost, `helios-python`
   the only pyo3 consumer. Full map in `ARCHITECTURE.md`.
4. **No speculative crates.** Only `helios-core` is a member now; the other 10 crates
   are created when their layer is implemented (architecture_scoping growth triggers).
   The full Atlas dependency set is declared in `workspace.dependencies` as SSOT.
5. **Foundation numeric representation.** `helios-core` constants/newtypes are `f64`
   at the definition boundary; the generic `Scalar` seam lands in `helios-math`.

## Delivered

- Workspace skeleton + config (`Cargo.toml`, `rust-toolchain.toml`,
  `.config/nextest.toml`, `.gitignore`).
- `helios-core` (0.0.1): typed errors, CODATA/ICRU constants with derivation tests,
  validating newtypes. 13 tests green; clippy `-D warnings` clean.
- Foundation artifacts: README, ARCHITECTURE, backlog, CHECKLIST, gap_audit,
  CHANGELOG.

## Metrics

| Metric | Value |
|--------|-------|
| Crates implemented | 3 / 11 (`helios-core`, `helios-math`, `helios-domain`) |
| Tests | 25 passed / 0 failed |
| Clippy warnings (production) | 0 |
| Test wall-clock | 0.23 s (well within 30 s budget) |

## Gaps opened

See `gap_audit.md` G-1..G-11 (physics, numerics, accuracy, integration, tooling).
Highest-risk: G-1 (no physics), G-3 (no accuracy oracles), G-5/G-11 (Atlas API
surfaces / geometry ownership).

## Correction (user directive)

Geometry primitives (`Aabb`/`Ray`) belong to **gaia**, not Helios. The first
`helios-math` cut defined them locally; they were removed as a downstream
duplication. gaia already owns `Aabb` and a validated `Ray`+`intersect_aabb`
(leto-migration branch). `helios-math` now exports only the `Scalar` seam and the
leto substrate; gaia geometry is consumed via H-003b once gaia's migration lands.

## Next increment

**H-004b:** `helios-domain` — verified `ritk-io` DICOM (CT/MVCT) load path into
`Volume`/`VoxelGrid` (pose from `ImageOrientationPatient`, HU rescale), plus
`CtVolume`/`MvctVolume` HU-semantic newtypes. Decomposed plan in `CHECKLIST.md`.
(ritk pulls burn+dicom — heavy build.)
