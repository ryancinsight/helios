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
| Crates implemented | 1 / 11 (`helios-core`) |
| Tests | 13 passed / 0 failed |
| Clippy warnings (production) | 0 |
| Cold build | 18.1 s |
| Test wall-clock | 0.19 s (well within 30 s budget) |

## Gaps opened

See `gap_audit.md` G-1..G-9 (physics, numerics, accuracy, integration, tooling).
Highest-risk: G-1 (no physics), G-3 (no accuracy oracles), G-5 (Atlas API surfaces
unverified against real usage).

## Next increment

**H-003:** `helios-math` — sealed `Scalar` seam over `hermes`/`leto` + geometry
primitives (`Vec3<T>`, affine patient/beam transforms, ray/AABB intersection) with
analytically-derived intersection tests. Decomposed plan in `CHECKLIST.md`.
