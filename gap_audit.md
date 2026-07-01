# Helios Gap Audit

Physics, numerics, accuracy, architecture, and integration gaps. Closed by
evidence, not silence. Each gap: ID, description, class, current evidence tier,
target closure.

## Open gaps

### Physics / numerics

- **G-1 (physics):** No radiation-interaction physics yet. Need photon linear
  attenuation coefficients (μ/ρ) from a citable source (NIST XCOM / ICRU) and
  electron transport model. *Evidence tier: none.* → H-011.
- **G-2 (numerics):** ~~No `Scalar` seam.~~ **CLOSED (H-003).** `helios-math`
  establishes `Scalar = eunomia::RealField` (the Atlas numeric SSOT) as the Helios
  compute seam; `Ray`/`Aabb` are generic over it and exercised at both `f64` and
  `f32` (native precision, no widen-narrow). `helios-core` constants remain `f64`
  literals by design and are converted by callers.
- **G-3 (accuracy):** No dose-engine or projector reference solutions yet. Gamma
  (3%/2mm), DVH, and MVCT image-quality oracles are unimplemented. Validation vs
  VoLO/TOPAS/GATE/EGSnrc pending. *Evidence tier: none.* → H-012, H-013, H-042.
- **G-4 (numerics):** Reduction-order sensitivity for future GPU vs CPU differential
  tests not yet characterized; epsilon bounds must be derived per reduction depth
  when the projector/dose kernels land. → H-012.

### Architecture / integration

- **G-5 (integration):** Atlas crate *APIs* partially exercised. **eunomia**
  (`RealField`/`FloatElement`/`NumericElement`) and **leto** (`Vector3`) now
  verified against real usage and building in-tree (H-003). `ritk-io` (DICOM/MVCT),
  `gaia` (MLC geometry), hephaestus, moirai, coeus, consus surfaces remain
  unverified; symbol existence must be confirmed via `cargo doc`/source before each
  first use (anti-hallucination). → H-004, H-005, H-010+.
- **G-10 (integration, upstream co-evolution):** leto's **default** features pull
  `mnemosyne` at a rev pinned to `themis ^0.8`, which conflicts with themis HEAD
  `0.9.17` — a version skew in the Atlas stack's transitive git graph. *Workaround
  applied:* Helios consumes leto with `default-features = false, features=["std"]`,
  deferring mnemosyne placement to the layer that needs it (themis/mnemosyne
  integration, later sprint). *Upstream item:* leto's pinned mnemosyne rev (or
  mnemosyne's themis bound) should be advanced to themis 0.9.x so the default
  feature set resolves. File against the leto/mnemosyne repos when that layer is
  built. *Evidence tier: reproduced (cargo resolution error), worked around.*
- **G-6 (build hygiene):** ~~Helios target-dir sharing.~~ **CLOSED.** Helios
  automatically routes its build through the shared `D:/atlas/target` via the
  inherited `repos/.cargo/config.toml` (`[build] target-dir`); Cargo discovers it by
  walking up from the package dir. Verified: `cargo doc` emitted to
  `D:/atlas/target/doc` and no per-`helios` `target/` exists. No action needed;
  backlog H-006 closed.
- **G-7 (toolchain):** `rust-toolchain.toml` pins `channel = "stable"` (currently
  1.95) but does not pin an exact version; MSRV floor declared as 1.85 in
  `Cargo.toml` (`rust-version`) but not yet CI-verified. → revisit at first CI.

### Testing / tooling

- **G-8 (coverage):** No `cargo-llvm-cov` run yet; >80% core-logic coverage target
  unmeasured. Applies from first `[minor]`. → measure after H-003.
- **G-9 (CI):** No CI pipeline wired (fmt/clippy/nextest/doc/audit/deny). Gates are
  run locally only. → file when the workspace has ≥2 crates.

## Closed gaps

- **G-6 (build hygiene):** Helios inherits the shared `D:/atlas/target` build dir
  via `repos/.cargo/config.toml`; no per-repo target. Verified this session.

## Residual risk register

- Atlas upstream APIs may drift (multi-repo co-evolution); Helios must pin commits
  in `Cargo.lock` and add cross-repo contract tests as it consumes each crate
  (G-5). Currently no lockfile committed for git deps because none are used yet.
- Physical constants (G-2) are CODATA-2018/ICRU-90 values verified by inter-constant
  derivation tests, not by an external authoritative fetch this session; values are
  standard and cross-checked, but a future audit should confirm against the live
  NIST database.
