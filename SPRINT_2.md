# Sprint 2 — Physics & GPU foundation

**Goal:** Begin the compute layers. First deterministic **photon-attenuation
physics** (the shared basis for MVCT forward projection and dose ray-tracing),
then the GPU foundation (hephaestus + moirai) and the first projection kernel.

**Phase:** 2 (Execution). Sprint 1's domain core (`helios-core`/`math`/`domain`)
is complete and gate-green.

## Sequencing note (unblocked vs blocked)

The goal acknowledges gaia's nalgebra/ndarray→leto/hephaestus migration as
concurrent, non-blocking context. Helios continues in parallel by ordering work so
geometry-independent physics lands first:

- **Unblocked now:** attenuation relations (H-011, done), NIST μ/ρ data (H-011b),
  CT calibration — pure material physics, no geometry.
- **Blocked on gaia geometry (G-11):** ray-marched line integral ∫μ dl (H-011c) and
  the GPU MVCT forward projector (H-012), which need `gaia::Ray` + voxel DDA once
  gaia's leto geometry reaches its default branch.
- **Blocked on heavy build:** ritk-io DICOM (H-004b) — sequenced as its own
  increment (ritk pulls burn wgpu+autodiff + dicom).

## Decisions (ADR-lite)

1. **Physics before geometry-coupled projection.** Attenuation coefficients,
   Beer–Lambert, and HU→density are implemented and analytically verified without
   any ray/voxel traversal, closing the highest-risk gap (G-1) immediately rather
   than waiting on gaia.
2. **Relations vs data separation.** `helios-physics` owns the *relations*
   (`μ=(μ/ρ)·ρ`, Beer–Lambert, calibration); concrete NIST XCOM μ/ρ *tables* are
   data loaded by H-011b. Tests assert the relations (analytical oracles), never a
   memorized cross-section digit — no fabricated reference values.

## Delivered (this increment)

- `helios-physics` (0.0.1): `LinearAttenuation`/`MassAttenuation` validated
  newtypes, Beer–Lambert `transmission`, `half_value_layer`, `to_linear`
  (`μ=(μ/ρ)·ρ`), and `relative_electron_density_from_hu`/`mass_density_from_hu`.
  9 analytical tests (`T(HVL)=½`, `T(0)=1`, density scaling, HU reference points,
  f32 genericity).

## Metrics

| Metric | Value |
|--------|-------|
| Crates implemented | 5 / 11 (`core`, `math`, `domain`, `physics`, `solver`) |
| Tests | 44 passed / 0 failed |
| Clippy warnings (production) | 0 |
| Test wall-clock | 0.45 s |

Also delivered: `helios-physics::projection` (line-integral reduction) and
`helios-solver::attenuation_map` (deterministic HU→μ engine — the first Sprint-2
compute kernel, CPU reference).

## GPU backend status (H-010, blocked — G-12)

Evidence-based finding: `hephaestus-wgpu` consumes the leto/mnemosyne/themis cluster
(path deps + `mnemosyne-memory` + pinned themis rev) — the same graph that failed
resolution in G-10 — plus a heavy `wgpu` build. Consuming it as a git dep would not
resolve cleanly against the current stack, which is mid-migration to a consistent
leto/hephaestus foundation (the migration the goal flags for gaia). Decision: author
every engine CPU-first; the GPU kernel is a differential drop-in against
`attenuation_map` once the stack converges. The `hephaestus_core::ComputeDevice`
seam and `hephaestus-wgpu` op surface are already scoped.

## Next increment

**H-011b:** energy-indexed NIST XCOM μ/ρ tables (water/air) + log-log interpolation
and material/mixture lookup, feeding `MassAttenuation` (values sourced/verified from
NIST). Unblocked. Decomposed plan in `CHECKLIST.md`.
