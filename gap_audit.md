# Helios Gap Audit

Physics, numerics, accuracy, architecture, and integration gaps. Closed by
evidence, not silence. Each gap: ID, description, class, current evidence tier,
target closure.

## Inverse-planning dose objective metrics (H-099, 2026-08-05)

The live planning audit found one remaining public Aequitas gap spanning the
autodiff-only `helios-planning` objective layer and the shared
`helios-analysis::Dvh` gEUD entry points: `DvhPenalty` exposed clinical dose
floors and ceilings as `&[f64]`, `EudPenalty` exposed its gEUD reference dose
as `f64`, and the public gEUD volume-effect parameters were untyped
dimensionless scalars. These were physical planning inputs, not dense tensor
storage.

H-099 closes the gap. DVH bands and gEUD references now use Aequitas
`AbsorbedDose<f64>`, and the gEUD parameter uses `Dimensionless<f64>`. Base
values are extracted only at the Coeus or Asclepius formula boundaries. Beamlet
weights, penalty coefficients, response slopes, and dense dose-influence
entries remain scalar because their units are optimization-model coefficients
rather than fixed SI metrics.

The planning law is real-valued under Eunomia. It has no phasor boundary, so no
imaginary dose unit or complex physical wrapper is introduced. The focused
all-feature package check passes; the value-semantic autodiff suite retains its
independent gEUD, finite-difference, hinge, and optimizer coverage. The
committed lock now includes the direct planning-to-Aequitas edge with Git
source identity, and the exact clean-source locked package check passes. Hosted
run `31011688127` passes the Rust, Python, dependency, and phase-replicated
benchmark gates at exact head `c00d270`; the classifier reports 0 regressions
and 0 replication-universe mismatches. See [`ADR 0017`](docs/adr/0017-planning-dose-quantities.md).

## Eunomia complex compatibility refresh (2026-07-28)

The live Helios source contains no `Complex`, `Complex32`, `Complex64`, or
imaginary-valued public physical contract. Its Aequitas migration therefore
needs no complex-unit extension or consumer wrapper. Eunomia's complex support
is relevant to other numerical consumers, but no Helios metric crosses that
boundary.

## Radon imaging geometry (H-097, 2026-08-02)

The live imaging audit found a remaining public Aequitas gap in
`helios-imaging`: `Sinogram` stored projection angles and detector offsets as
raw `T`, while Radon and SIRT accepted raw source distance and ray-march step.
The documentation named radians and millimetres, but the type system did not
enforce those dimensions.

H-097 carries `Angle<T>` and `Length<T>` through `Sinogram`, Radon, FBP, SIRT,
the quantum-noise path, all affected examples, the analysis validation example,
and the end-to-end workflow. Scalar extraction is confined to radians for
trigonometry, millimetres for ray/grid geometry, and centimetres for the FBP
ramp filter. Dense voxel-grid coordinates and line-integral readings remain
scalar kernel storage.

This imaging law is real-valued. It has no phasor or Fourier-valued physical
field, so Eunomia complex values do not require an imaginary length unit or a
consumer wrapper. Complex compatibility remains the existing Eunomia rule:
complex components share one real physical unit when a numerical domain
actually needs a phasor.

The earlier provider rev-pin removal left the committed lockfile without the
current git source identities. H-097 refreshes that derived lock state so the
standalone `--locked` gate is executable again.

Evidence: `cargo check --offline --locked -p helios-imaging --all-targets`,
warning-denied Clippy, `cargo test --doc`, Rustdoc, and the analysis
`validation_regression` example check pass. Nextest run
`26905c71-03f1-447a-be77-df0c84278c3c` passes all 33 `helios-imaging` tests,
including the analytical disk oracle and f32/f64 reconstruction paths.

The downstream `cargo nextest -p helios-analysis -p helios-simulation` gate
does not reach test execution: `mnemosyne-memory` exits from `rustc` with code
1 and no diagnostic. The simulation example check reaches the provider graph
but is rejected by Moirai's existing `missing_docs` errors for public fields in
`moirai-executor/src/metrics/mod.rs` and `moirai-executor/src/task/mod.rs`.
These are external provider integration blockers, not failures of the typed
Radon contract; they remain open in the provider-owned audit path.

## Historical benchmark and provider-graph gate (H-098, 2026-08-04)

The H-097 hosted benchmark job reached candidate compilation but failed while
loading the historical baseline workspace: its committed `Cargo.toml` still
declares `../moirai/moirai` and `../moirai/moirai-parallel`, and the clean
runner did not contain those sibling checkouts. The failure occurred before a
benchmark binary or metric was executed (`30913557127`).

The corrected workflow materializes the candidate and historical baseline path
dependencies through the pinned Atlas checkout tool at workspace `.`. This is
required because the baseline's `../moirai/moirai` and
`../moirai/moirai-parallel` paths resolve relative to its checked-out workspace;
using `..` placed the repositories one level too high and caused the checkout
action to inspect a non-repository directory. Baseline metadata and both
benchmark phases require `--locked`, and the job no longer marks this failure
class successful with `continue-on-error`.

The first replacement run also showed that an additional `coeus/Cargo.toml`
checkout invocation was invalid: Helios consumes Coeus as a Git dependency and
the job contains no `coeus/` checkout, so the action failed before the Atlas
benchmark gate. The invocation is removed; Coeus remains resolved through the
candidate's locked Git graph rather than a fabricated path checkout.

The first clean lock regeneration then exposed a separate provider defect:
Gaia and Asclepius still required Eunomia `^0.7.0` while current Aequitas and
Leto require Eunomia `0.8.0`. Gaia PR `#21` and Asclepius PR `#6` corrected
those manifests and are merged at `683565e` and `3463a70`. Helios now uses
their normal Git dependencies, and `Cargo.lock` was regenerated outside
Atlas's local path overlay so its first-party entries carry real Git source
identities. No temporary pin or compatibility shim remains.

This is a verification-infrastructure and dependency-graph correction, not a
new physical metric or unit. It does not change the benchmark source,
workload, counterbalancing, classifier, or Eunomia boundary. The replacement
hosted matrix is the closure oracle; no Aequitas gap is inferred from the
pre-execution failures. Helios's public physical contracts remain real-valued;
no imaginary or complex unit is required. Local locked metadata,
warning-denied workspace Clippy, 285/285 configured Nextest tests, doctests,
and warning-denied Rustdoc pass against the clean graph; the post-merge-graph
Nextest run is `0eff8ef7-ca1e-4e50-9a8e-c9e070e671cc`. Hosted run
`31011688127` passes all required jobs with 0 benchmark regressions and 0
replication-universe mismatches. H-098 is closed.

## Aequitas metric audit refresh (2026-07-27)

The current `main` branch had stale closed-ledger claims: the typed helical
delivery and collimation implementations existed on a non-ancestor branch but
were absent from the live source. The omission reintroduced raw physical
scalars in public contracts. H-093 restores the existing provider-first
implementation and migrates the current branch's callers without a scalar
compatibility path.

- `helios-domain::HelicalDelivery` stores and returns Aequitas `Length`,
  `Dimensionless`, `Time`, `Angle`, and `Velocity`; `FieldAperture` stores its
  penumbra as `Length`.
- `helios-simulation::HelicalProjection`, `DeliveryFrame`, acquisition
  geometry inputs, and delivery/dose/portal angle boundaries now preserve the
  typed values through public APIs. Conversion to millimetres/radians occurs
  only at the Gaia/scalar trigonometry kernels.
- Focused `cargo check`, warning-denied Clippy, configured Nextest, doctest,
  and direct rustfmt checks completed offline. Cargo emitted only the existing
  unused local-patch warnings from the dirty shared provider manifest; those
  manifest and lockfile changes are outside H-093 and remain uncommitted.

The remaining audit boundary is unchanged: scalar dense voxel storage,
optical-depth/gamma values, fractions, indices, and probabilities are not
dimensional public metrics. H-093 is closed by the source/type audit and the
focused value-semantic gates; a full locked workspace gate still depends on
the shared provider graph represented by the pre-existing dirty manifest and
lockfile.

## Scatter-kernel physical inputs (H-095, 2026-07-31)

The audit found a public Aequitas gap in `helios-solver`: exponential,
symmetric, forward-peaked, poly-energetic, and oriented scatter APIs accepted
physical ranges and sampling pitches as raw `T`, while
`SpectralComponent` also exposed a raw relative weight. `helios-simulation`
repeated the untyped ranges in `CollapsedCone` constructors and stored the
oriented sampling step as a scalar.

H-095 closes the gap by carrying `Length<T>` through all public range and
sampling-pitch boundaries and `Dimensionless<T>` through spectral weights.
Scalar conversion remains only at the centimetre exponential formula and
millimetre voxel-mesh boundaries; dense taps and dose fields remain scalar
storage. The public unit-suffixed names were removed rather than retained as
forwarding shims. Eunomia compatibility is real-only: these transport kernels
have no phasor, Fourier, or imaginary-valued physical contract, so no complex
quantity extension is required.

Evidence: targeted `cargo check -p helios-solver -p helios-simulation
--all-targets --offline` passes after migrating tests, examples, and the
benchmark. Focused Nextest, warning-denied Clippy, doctest, Rustdoc, and the
full locked-workspace gate remain closure checks. The latter also depends on
the pre-existing dirty shared provider manifest/lockfile and peer workspace
changes, which are outside this item.

## Live GPU attenuation refresh (2026-07-27)

`HELIOS-AEQ-MET-05` closed a live public boundary missed by the earlier audit:
`helios-gpu::GpuAttenuationMapper::new` accepted mass attenuation and water
density as raw `f32` values even though the CPU attenuation path already used
typed Aequitas/provider quantities. The constructor now accepts Aequitas
`AreaPerMass<f32>` and `MassDensity<f32>`; conversion to the cm-based GPU
kernel units occurs once at the explicit formula boundary. Tests and the GPU
example use the typed constructors, with no scalar compatibility facade.

`helios-gpu` check with tests/examples, Nextest (10/10), warning-denied
Clippy, doctests, Rustdoc, rustfmt, and diff checks pass against the shared
offline graph. Existing unused local-patch warnings and the dirty peer
lockfile remain outside this metric slice.

## Aequitas metric gap audit (2026-07-23)

The dose field itself remains `helios_domain::Volume<T>` storage. This audit
targets values that cross a public analysis, delivery, or geometry API. It does
not count dense field storage, gamma values, volume fractions, homogeneity
indices, or response probabilities as dimensional gaps.

### Existing coverage

Helios already uses Aequitas `EnergyPerArea` for portal fluence, `AbsorbedDose`
for deposition totals and DVH samples, `Length`/`ReciprocalLength` for geometry
and attenuation, and `AreaPerMass` for mass attenuation. `EnergyMeV` and
`VoxelSpacingMm` are validated consumer newtypes backed by Aequitas quantities.

### Open implementation ledger

| ID | Evidence | Remaining implementation | Owner | Status / acceptance oracle |
|---|---|---|---|---|
| `HELIOS-AEQ-MET-01` | `helios-analysis/src/dvh.rs` stores `Vec<AbsorbedDose<T>>`, but `min`, `max`, `mean`, `dose_at_volume_fraction`, and `generalized_eud` return raw `T`; dose criteria and TD50/TCD50 parameters also enter as `T`. | Return `AbsorbedDose<T>` for dose-valued results and parameters; keep Vx, HI, TCP, and NTCP dimensionless or probability-typed. | Helios | **RESOLVED.** `Dvh` extrema, mean, Dx, gEUD, and TCP/NTCP dose parameters now use `AbsorbedDose<T>`; nearest-rank, masked, NaN, Asclepius-law, and end-to-end PTV/OAR value semantics remain covered. |
| `HELIOS-AEQ-MET-02` | `helios-analysis/src/gamma.rs` accepted `dta_mm`, normalization dose, low-dose cutoff, and search radius as raw `T`; gamma volume and pass rate are dimensionless. | Type distances as `Length`, dose thresholds as `AbsorbedDose`, and keep the result storage scalar/dimensionless. | Helios | **RESOLVED.** `gamma_index_3d`, `gamma_index_3d_local`, and `gamma_pass_rate` now type physical criteria with Aequitas while retaining the Low gamma kernel, local/global normalization, grid checks, scalar gamma field, and scalar pass rate. Focused value-semantic gamma tests and all in-tree callers migrate; ADR 0007 records the breaking boundary. |
| `HELIOS-AEQ-MET-03` | `helios-simulation/src/delivery.rs` stored leaf fluence as `T`; `total_delivered_fluence` returned `T`. `portal.rs` constructed `EnergyPerArea` internally and converted it back. `dose_accumulation.rs` accepted `*_mm` geometry and sampling values as `T`. | Carry fluence as `EnergyPerArea` and geometry distances as `Length` through delivery, portal, and dose accumulation. | Helios | **RESOLVED.** `DeliveryFrame`, collimation, portal transmission, total fluence, and dose geometry now use Aequitas quantities. Typed values convert once to the existing millimetre ray/voxel kernel; closed-leaf zero, Beer–Lambert darkening, fluence linearity, geometry-limit, f32, example, and end-to-end regressions pass. ADR 0008 records the breaking boundary. |
| `HELIOS-AEQ-MET-04` | `helios-analysis/src/image_quality.rs` returns raw intensity/RMSE values, while the same analysis can operate on dose volumes. | Decide the semantic input at the analysis boundary: retain raw image intensity for MVCT, but return `AbsorbedDose` statistics when the API contract is dose-specific. | Helios | **RESOLVED.** Shared ROI/RMSE value kernels now back raw MVCT `roi_statistics`/`volume_rmse` and typed-dose `dose_roi_statistics`/`dose_volume_rmse`; the clinical validation example uses typed dose means/stddev and converts only for dimensionless contrast/CNR. Value tests cover f64/f32 and Gray outputs; ADR 0009 records the partition. |
| `HELIOS-AEQ-MET-05` | `helios-gpu/src/attenuation.rs` accepted mass attenuation `mu_over_rho_cm2_g` and reference water density `water_density_g_cm3` as raw `f32` at `GpuAttenuationMapper::new`, despite the CPU attenuation boundary using Aequitas/provider quantities. | Accept typed Aequitas `AreaPerMass<f32>` and `MassDensity<f32>`; extract cm²/g and g/cm³ only at the GPU formula boundary and migrate tests/examples without a scalar facade. | Helios, Aequitas | **RESOLVED in this increment.** `GpuAttenuationMapper::new` carries typed coefficient/density inputs and preserves the fused clamp law. Check with tests/examples, Nextest 10/10, warning-denied Clippy, doctests, Rustdoc, rustfmt, and diff checks pass. See [ADR 0013](docs/adr/0013-gpu-attenuation-quantities.md). |

| `HELIOS-AEQ-MET-06` | `helios-imaging::Sinogram`, `parallel_beam_radon`, and `sirt_reconstruction` accepted projection angles, detector offsets, source distance, and ray step as raw `T` despite their radian/millimetre contract. | Carry Aequitas `Angle<T>` and `Length<T>` through the public imaging boundary; extract only at trigonometry, mesh, and filter formula boundaries. | Helios, Aequitas | **RESOLVED in H-097.** Radon, FBP, SIRT, noise, examples, analysis validation, and the end-to-end workflow use typed geometry; the standalone locked gate and focused value-semantic suites are the closure evidence. See [ADR 0016](docs/adr/0016-radon-geometry-quantities.md). |
| `HELIOS-AEQ-MET-09` | `helios-physics` Compton APIs accepted raw scalars documented as MeV for photon energy, so the public boundary carried no energy dimension. | Accept Aequitas `Energy<T>` for Klein–Nishina, energy-transfer, and Compton mass-coefficient APIs; convert to MeV only at the dimensionless kernel boundary. | Helios | **RESOLVED in H-101.** The Rust boundary uses typed Aequitas energy, the Python adapter retains validated MeV input only at FFI, and the 1 MeV/1,000,000 eV value-equivalence oracle passes. Current Atlas metadata resolves `coeus-autograd` through `repos/coeus/crates/coeus-autograd`; full warning-denied Clippy, 286/286 Nextest, doctest, Rustdoc, format, and package Python gates pass. See [ADR 0010](docs/adr/0010-compton-energy-quantity.md). |

### Explicit non-gaps and constraints

- TERMA/dose arrays remain scalar storage until a field-descriptor contract can
  carry dimensions without changing dense-kernel representation.
- Beam angles, gamma values, fractions, homogeneity indices, and TCP/NTCP are
  dimensionless; they must not be wrapped as length or dose merely because they
  are reported beside physical quantities.
- `HELIOS-AEQ-MET-04` is closed. Raw MVCT intensity remains scalar, while
  dose-specific ROI and RMSE results carry `AbsorbedDose`; no Aequitas provider
  extension was required. Future public signature changes must update their
  examples, Python surface, and focused tests in the same change.

## H-003d oriented-grid provider convergence (closed)

- Leto 0.38 now owns checked conversion from world-space rotation columns to a
  `UnitQuaternion`; Helios uses that provider contract to restore the
  local-index-to-world `Isometry3` grid pose. CPU ray clipping and terma
  deposition run in scaled-index space and retain the world-space millimetre
  parameter; the HDF5 boundary records three validated rotation columns. The
  grid core does not duplicate DICOM tags or matrix-to-quaternion logic.
- The current Hephaestus `FieldGeometry` has no rigid pose. `GpuProjector`
  therefore returns a typed dispatch error for a non-identity rotation before
  data upload (correctness evidence: type/contract plus a value-semantic test),
  rather than silently omitting orientation. A pose-bearing GPU field geometry
  remains an upstream Hephaestus capability gap.
- Evidence tier: type-level rigid pose and checked Leto basis construction;
  analytical/differential validation through the 104/104 focused nextest run
  (including oriented Beer–Lambert and HDF5 pose round trips, live GPU checks);
  warning-denied Clippy, doctest/rustdoc, workspace example build, workspace
  format check, and four 196/196-package SemVer checks are clean.
- H-004d remains externally sequenced: RITK's public DICOM tags currently omit
  `ImageOrientationPatient`, and both permitted RITK worktree lanes carry
  active peer migrations. Helios will consume the named provider tag once that
  owner lane is available.

## Open gaps

### H-088 — deterministic book-figure SSOT gate (implemented, PR #32 merged)

- `xtask` now owns `FIGURE_SPECS`, deterministic SHA-256 manifest generation,
  and the `check-figures` command. The command validates the seven committed
  SVGs, scans `docs/book/SUMMARY.md` and `docs/book/README.md`, and fails on
  either an unlisted asset or a docs/spec mismatch.
- The workflow runs the gate after Rustdoc. The Python and benchmark lanes
  allow Cargo to refresh the lock after Atlas path-dependency materialization;
  the Rust workspace remains locked. The previous CI failure was the absent
  subcommand, not a figure mismatch.
- Local evidence: `cargo check -p xtask --offline`, formatter check, `mdbook
  build docs/book`, and `cargo run -p xtask --offline --locked -- check-figures`
  pass with `SSOT_IN_SYNC` and 7/7 references. Hosted PR #32 run
  `30070400660` passes build, Rust workspace, Python bindings, and the
  replicated benchmark gate; the PR merges as `02d7a775`.

### H-087 — portal fluence quantity boundary (implemented, PR #32 merged)

- `helios-simulation::frame_portal_fluence` now carries the transmitted portal
  fluence as Aequitas `EnergyPerArea<T>` through Hyperion's dimensionless
  transmission product before converting at the established scalar frame API.
  The direct Aequitas pin is `3ae0b6b`; implementation commit `b2a9ebe`.
- This closes the remaining quantity-conversion seam in the portal workflow;
  the dense fluence frame remains representation storage, not a second metric
  owner. The direct provider pin is now the merged Aequitas revision
  `e0fc5f3`. Existing full-transmission, Beer–Lambert, closed-leaf, f32, and
  invalid-optical-depth regressions remain the behavioral oracle.
- PR #32 merges as `02d7a775` from implementation head `31147f0` and PM
  follow-up `5832ffa`. The hosted build, Rust workspace, Python, and
  replicated benchmark checks pass. A local focused `helios-simulation`
  Nextest attempt could not start because active peer CFDrs/Leto builds held
  the shared Atlas lock; no local package-gate result is claimed.

### G-29 — DICOM charset dependency (externally blocked)

- `dicom-encoding` 0.10.0 declares `encoding` 0.2.33 unconditionally and uses
  it for the DICOM Specific Character Set codecs. The current release exposes
  no feature that can remove the dependency. A 2026-07-20 registry and resolved
  manifest audit confirms 0.10.0 remains the latest release and still carries
  the unconditional edge.
- RUSTSEC-2021-0153 reports maintenance status, not a known vulnerability.
  CI quarantines only that advisory ID while continuing to deny every other
  warning and vulnerability. Reopen H-073 when `dicom-rs` publishes a release
  backed by a maintained charset implementation; a consumer fork or reduced
  character-set implementation would duplicate provider ownership.

### Recently closed

- **G-30 — RESOLVED (H-083).** A stale unclaimed mdBook expansion described
  numerous types and methods absent from the current Helios source, while four
  older appendix/numeric pages contained control characters created by escaped
  Markdown delimiters. The expansion was not published. The recovered book now
  names APIs verified against crate re-exports, delegates release history and
  architecture to their root SSOT documents, contains no control characters,
  resolves every relative Markdown link, and builds without warnings under
  mdBook 0.5.4. The legacy scanner independently proves that its emptied
  migration allowlist matches zero current legacy dependency or source surfaces.

- **G-27 — RESOLVED (H-071).** The copied same-run classifier and bare native
  test invocation are deleted. Implementation head `44fb2768d` uses Atlas gate
  `9bfb722`, holds the candidate Criterion harness constant, measures
  phase-reversed ABBA and BAAB replications, runs native tests through the
  committed Nextest budget, and passes Rust, Python, and benchmark jobs in
  hosted run `29784712768`.

- **G-28 — RESOLVED (H-072).** The isolated Python binding crate previously
  resolved PyO3 0.23.5, which is affected by RUSTSEC-2025-0020 and
  RUSTSEC-2026-0177. PyO3 0.29.0 closes both vulnerabilities; the thin boundary
  uses `Python::detach` around the existing Rust planning call and retains no
  domain logic. The built extension is covered by the value-semantic Python
  contract suite. G-29 records the sole exact unmaintained dependency
  quarantine.

- **G-26 — RESOLVED (H-068).** `EnergyMeV` and `VoxelSpacingMm` previously
  stored dimensionless `f64` values despite representing physical quantities.
  Their validated newtypes now store Aequitas `Energy<f64>` and `Length<f64>`,
  preserve MeV and millimetre at the public boundary, and retain their
  zero-overhead scalar layout through compile-time size/alignment assertions.
  Round-trip properties use a bound derived from four machine-epsilon
  roundings. `HounsfieldUnit` remains Helios-owned because it is a calibrated
  non-SI scale. Evidence: warning-denied all-target Clippy, 17/17 configured
  Nextest tests, doctests, and warning-clean rustdoc for `helios-core`.

- **G-25 — RESOLVED (H-067).** The stale local lock changed only the
  `apollo-fft` version field, which was not a complete Cargo resolution and
  failed the warning-denied locked gate. Regenerating the Apollo package
  closure selects `apollo-fft` 0.25.0, Eunomia 0.4.0, Leto 0.38.2, and
  Hephaestus 0.17.0, removes Eunomia's obsolete `num-traits` edge, and removes
  Hephaestus WGPU's `num-complex` edge plus the package itself. The root
  manifest already follows the Apollo default branch; no source, manifest,
  compatibility wrapper, or fallback change is required. Locked metadata and
  format pass; warning-denied all-target/all-feature workspace Clippy passes;
  configured Nextest is 272/272; all ten Rust library doctest targets pass
  with zero examples; workspace rustdoc is warning-clean. Evidence tier:
  compiler-checked dependency resolution, warning-denied diagnostics, and
  value-semantic workspace regression execution.

- **G-24 — RESOLVED (H-066).** The workspace declared `num-traits` directly
  even though no Helios manifest or source consumed it. The direct declaration
  is removed. `cargo check --workspace --locked` passes, and
  `cargo tree -i num-traits --locked --edges normal` shows only transitive
  provider paths through Eunomia, Gaia, Half, WGPU, and their dependencies.
  Evidence tier: compiler-checked dependency resolution plus an inverse
  dependency-tree ownership audit.

- **G-23 — RESOLVED (H-005 reconciliation).** The foundation roadmap still
  listed a binary-MLC plus collimator/jaw model as todo after H-020b delivered
  `LeafOpenTimeSinogram`/`MlcModel` and H-020k delivered gaia-backed
  `FieldAperture` plus `collimate_frames`. The board marks duplicate H-022
  done and the README states the delivered ownership without creating a
  duplicate implementation track.
  The same reconciliation check found and restored the exact workspace formatter
  output in `helios-math::lib` and the solver deposition/projector tests.
  Evidence tier: source/API and board reconciliation plus formatter verification.

- **G-19 — RESOLVED (H-062).** Repeated `Dvh::volume_fraction_at_dose` queries
  previously scanned the complete sorted sample for every threshold, which made
  a plan with `q` DVH queries O(q·n). The query now uses a zero-allocation
  `partition_point` lower bound over the existing sorted slice, reducing the
  finite/infinite path to O(q·log n). A `contains_nan` marker preserves the
  previous `>=` filter semantics for unordered samples, and a NaN threshold
  returns zero as before. The fixed 64³/1,024-query Criterion workload measured
  30.090 ms [29.717, 30.472] for the scan reference and 29.229 μs [28.426,
  30.023] for production after the change; the paired median ratio is 1,029×.
  Value-semantic focused nextest passes 34/34. Evidence tier: empirical
  Criterion comparison plus value-semantic/differential boundary tests; see
  `validation_reports/2026-07-15-dvh-query-optimization.md`.

- **G-20 — RESOLVED (H-063).** Helios's direct DICOM dependency declared 0.8
  while the local `ritk-dicom` provider supplied `dicom-object` 0.10. The
  resulting duplicate `dicom_core::Tag` types caused four E0308 errors in
  `helios-domain/src/dicom.rs` during the workspace example check. Helios now
  declares DICOM 0.10, and the lockfile resolves one DICOM 0.10.0 graph across
  Helios and `ritk-dicom`. Locked workspace examples, all-target all-feature
  Clippy, 261/261 workspace nextest tests, doctests, and rustdoc pass. Evidence
  tier: compile-time dependency/type verification plus value-semantic tests.

- **G-21 — RESOLVED (H-064).** Helios previously used `dicom::core` and
  `dicom::object` directly for typed attributes and synthetic test fixtures,
  despite `ritk-dicom` owning the DICOM boundary. Helios now consumes only the
  `ritk-dicom` public tags, attribute-read trait, parser, transfer-syntax, and
  decoder contracts. Production and test scans contain no direct `dicom::`
  imports; the focused provider-backed domain suite remains 41/41.

- **G-22 — RESOLVED (H-065).** Helios's lockfile still selected Moirai 0.2.0
  after the upstream 0.3.0 release retired its unowned NUMA iterator and
  benchmark. The regenerated lockfile selects 0.3.0 for every Moirai package;
  `cargo check --workspace --examples --all-features` compiles the complete
  example graph. Evidence tier: compiler-checked dependency resolution and
  example compilation.
  Evidence tier: dependency/identifier scan plus value-semantic nextest.

- **G-14 — RESOLVED (H-003c).** The concurrent leto geometry rewrite settled: leto
  and gaia now build against the new `leto::geometry` API (Vector3/Point3 with
  `.x/.y/.z` fields; `Isometry3` reduced to `{rotation, translation}`). Helios was
  adapted: `helios-math` re-exports the new leto types (`Point2/Point3/Vector3/
  UnitVector3`) + gaia `Aabb/Ray`; `VoxelGrid` simplified to **axis-aligned**
  (origin + spacing), dropping the now-incomplete `Isometry3` pose (oriented grids
  are a follow-up pending a rigid-transform primitive with `transform`/`inverse`);
  the projector's pose-rotation check was removed. **Full workspace builds; 97 tests
  pass** (all crates incl. live GPU), clippy `-D warnings` clean, fmt clean. The
  previously-blocked dose kernel-superposition engine (H-013b) is now built and
  verified. *Evidence tier: verified — whole-workspace build + 97 tests green.*

### (historical) BLOCKER — concurrent foundation refactor

- **G-14 (integration, BLOCKED — concurrent leto geometry relocation):** Mid-session
  the shared **leto** submodule advanced (peer/concurrent work) and its `geometry`
  module (`Vector3`/`Point3`/`Isometry3`/`UnitQuaternion`/…) is no longer present at
  leto's current HEAD (`git ls-files crates/leto/src/geometry` is empty; earlier this
  session gaia built 927 tests against `leto::geometry`). The types are not yet
  relocated to a discoverable home, so **gaia fails to compile** (86 errors,
  `unresolved import leto::geometry`), and every Helios crate that depends on
  `helios-math` (which re-exports gaia geometry) transitively fails to build —
  including at the last green commit `2ce36787` (the foundation shifted under it).
  *Interpretation:* the Atlas stack is mid-migration to **gaia-native geometry**
  (geometry moving out of leto into gaia), the end-state the earlier feedback
  intended. *Action (discipline):* do NOT fix leto/gaia's in-flight relocation
  (peer's active, cross-stack work; unknown target), do NOT revert the shared
  submodule, do NOT commit broken/unverified Helios code. **Deferred item H-013b
  (dose kernel superposition, `dose_convolution_x` + `exponential_deposition_kernel`)
  is written in `crates/helios-solver/src/dose.rs` with exact analytical oracles
  (delta-kernel identity, normalized-kernel interior conservation, physical build-up)
  but is UNVERIFIED — it cannot build until the geometry foundation settles.**
  *Unblock:* when gaia's native geometry lands, update `helios-math` to re-export all
  geometry from gaia (H-003c), then verify + commit H-013b. *Evidence tier:
  reproduced (leto HEAD has no geometry; gaia 86-error build failure).*
  *Mitigation (H-055):* `helios-math`'s geometry vocabulary is now behind a default
  `geometry` feature and `helios-physics` consumes it with `default-features=false`,
  so `helios-core`, `helios-math` (scalar seam) and `helios-physics` **build/test
  independently** of the churning geometry stack (`cargo nextest run -p helios-core
  -p helios-physics`). Only geometry-dependent crates (`helios-domain`/`-solver`,
  whole-workspace `cargo build`) remain blocked until the foundation settles.

- **G-15 (imaging accuracy):** *Partially addressed (H-030, H-033).* MVCT
  reconstruction (parallel-beam FBP) validated by a forward→reconstruct round-trip on a
  disk phantom, now *quantified* with `helios-analysis::image_quality` metrics
  (interior-ROI accuracy within 15% of μ₀, background suppression, disk/air contrast
  >0.85, CNR >1), and *quantum noise* (H-033b: `Var(τ')≈e^{τ}/N₀` validated vs
  analytical photon statistics; end-to-end noisy-recon noise scales with flux).
  **Remaining:** statistical reconstruction (OS-SEM/MLEM) and validation vs *published
  TomoTherapy MVCT data*; SIRT iterative reconstruction landed (H-030c, converges to its
  forward model, robust to noise/sparse-angle). The DICOM real-input path now ingests
  both a single slice
  (H-004b) and a full multi-slice **series** → 3-D HU `Volume` (H-004c:
  `load_ct_series` via `ritk-dicom`), so a real CT/MVCT study can drive the pipeline —
  clinical *dataset* validation still needs a licensed reference dataset. *Evidence tier: analytical/round-trip + synthetic-phantom
  metrics + real DICOM parse (synthetic round-trip through the ritk-dicom provider) — published-data
  comparison pending.*

### Physics / numerics

- **G-1 (physics):** *Closed through H-011 and H-011b.* Photon attenuation **relations**
  implemented and analytically verified in `helios-physics`: Beer–Lambert
  transmission, half-value layer, `μ = (μ/ρ)·ρ`, and first-order HU→density CT
  calibration (property/value-semantic tests: `T(HVL)=½`, `T(0)=1`,
  water/air/bone calibration points, f32 genericity). H-011b adds the selected
  NIST dry-air, liquid-water, and cortical-bone mass-attenuation knots over
  10 keV–20 MeV with an explicitly bounded log-linear interpolation contract.
  The five table value, boundary, interpolation, and invalid-domain tests pass
  on the final head. An electron-transport model remains a separate algorithm
  item, not a mass-attenuation table residual. *Evidence tier: analytical
  relations plus primary-source table-value tests.*
- **G-2 (numerics):** ~~No `Scalar` seam.~~ **CLOSED (H-003).** `helios-math`
  establishes `Scalar = eunomia::RealField` (the Atlas numeric SSOT) as the Helios
  compute seam and re-exports the leto linear-algebra substrate. `helios-core`
  constants remain `f64` literals by design and are converted by callers. The seam
  is exercised natively (`f32`/`f64`) by the first compute kernels as they land.
- **G-3 (accuracy):** *Partially closed (H-032).* The **validation machinery** now
  exists: `helios-analysis` implements the cumulative DVH (Dx/Vx/mean) and the 3D
  gamma index (Low, global normalization) + pass rate, with analytical oracles
  (identical→γ=0, criterion-scaled γ, uniform-DVH step, ramp quantiles). **Still
  open:** the dose-engine/projector *reference solutions* to validate (need
  H-013) and clinical comparison vs VoLO/TOPAS/GATE/EGSnrc (H-042). *Update:* the
  **MVCT forward projector** (`helios-solver::forward_project_ray`, H-011c) now
  produces line-integral projections (∫μ dl), analytically verified (uniform slab
  τ=μ·L, affine-field midpoint-exact). The dose engine's **primary-transport stage**
  (H-013a, `primary_fluence_parallel_x`) now produces the analytical exponential
  depth curve `Ψ₀·exp(−μx)`; the remaining stage is kernel superposition → dose
  (H-013b). Clinical comparison vs VoLO/TOPAS/GATE/EGSnrc (H-042) still pending real
  inputs.
  *Evidence tier: analytical (metrics + projector + primary transport verified) —
  full dose engine + clinical inputs pending.*
- **G-13 (numerics, units):** ~~Projector optical depth was 10× too large.~~
  **CLOSED (H-013a):** `forward_project_ray` mixed `μ` [cm⁻¹] with mm path length;
  now converts mm→cm for a dimensionless `τ`. Units contract documented (μ volumes
  in cm⁻¹, grid in mm). *Evidence tier: analytical (τ = μ·L_cm verified).*
- **G-4 (numerics):** Reduction-order sensitivity for future GPU vs CPU differential
  tests not yet characterized; epsilon bounds must be derived per reduction depth
  when the projector/dose kernels land. → H-012.

### Architecture / integration

- **G-5 (integration):** Atlas crate *APIs* partially exercised. **eunomia**
  (`RealField`/`FloatElement`/`NumericElement`) and **leto** (`Vector3`, `Point3`,
  `Isometry3`, `Translation3`, `UnitQuaternion`/`Unit`, `Array3` C-contiguous +
  `as_slice`) verified against real usage and building in-tree (H-003, H-004).
  `ritk-io` (DICOM/MVCT), `gaia` (MLC geometry), moirai, coeus, consus surfaces
  remain unverified; **hephaestus** `ComputeDevice` seam (GAT `Buffer<T: Pod>`,
  `alloc_zeroed`/`upload`/`download`/`write_buffer` with `themis::PlacementHint`)
  read and scoped for H-010 (not yet built — heavy wgpu compile + GPU-device
  availability are the gating risks). Symbol existence must be confirmed via
  `cargo doc`/source before each first use (anti-hallucination). → H-004b, H-005,
  H-010+.
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
- **G-7 (toolchain):** Helios 0.1.0 declares Rust 1.95, matching the merged
  provider graph. `rust-toolchain.toml` remains `stable` rather than an exact
  channel pin; the configured Rust 1.95 package gates are the current evidence.

- **G-11 (integration, geometry ownership):** Geometry primitives (`Aabb`, `Ray`,
  intersection, meshes, CSG) are owned by **gaia**, not Helios. gaia already has
  `Aabb<T: Scalar>` (over `leto::Point3`) and a validated-`UnitVector3` `Ray` with
  `intersect_aabb`. *Update (this session):* gaia's leto/eunomia migration is now
  **finalized and green** — gaia builds across all targets, **927 tests pass**,
  doctests pass, fmt clean; `Ray`/`Aabb` are committed and re-exported from `gaia`'s
  crate root (commits `b058eb0`, `ecd4060`). The source blocker is **resolved**.
  *Action taken earlier:* removed the duplicate `Ray`/`Aabb` from `helios-math`
  (upstream ownership). **Remaining (consumption wiring, H-003b):** the migration
  lives on gaia's `refactor/migrate-to-leto-geometry` branch, not yet merged to
  gaia's default branch — merging is a `refactor!` breaking change that also affects
  kwavers (co-evolution). **Update: consumption wired (H-050).** Helios now `[patch]`-
  redirects `leto`/`eunomia`/`gaia` git sources to the local synchronized checkout
  (one consistent source) and `helios-math` re-exports `gaia::{Aabb, Ray}`; a bridge
  test (gaia `Ray` ∩ gaia `Aabb` through Helios) passes. **G-11 is effectively
  closed** for local development — the projector (H-011c) is unblocked. Remaining
  release step: merge gaia migration to its default branch + update kwavers, then
  drop the patch. *Evidence tier: verified — Helios builds + 60 tests with local
  gaia geometry.*

- **G-12 (integration, GPU backend blocked):** `helios-gpu` on `hephaestus-wgpu` is
  blocked on the Atlas stack's leto/hephaestus dependency convergence — the same
  migration the goal flags ("gaia will need to move to leto/hephaestus"). Evidence:
  hephaestus's workspace consumes `leto`/`mnemosyne`/`themis` via **local path deps**
  with the `mnemosyne-memory` feature and a pinned `themis` rev, i.e. the same
  leto→mnemosyne→themis cluster that failed resolution in G-10, now compounded by a
  heavy `wgpu` build. Consuming `hephaestus-wgpu` as a git dep would not resolve
  cleanly against the current stack. *Decision:* do not force the GPU backend now;
  author every engine as a CPU reference first (`helios-solver`) so the GPU path
  (H-010) is a differential drop-in once the stack stabilizes. The
  `hephaestus_core::ComputeDevice` seam and `hephaestus-wgpu` op surface
  (`WgpuDevice::try_default`, `unary/scalar_elementwise_strided`, `reduction`) are
  already scoped for that increment.
  *Update (this session):* hephaestus is **verified green locally** — the workspace
  builds, `hephaestus-core` (21 tests) and `hephaestus-wgpu` (109 tests) pass, fmt
  clean, 0 code clippy warnings. Crucially the **wgpu GPU contract tests pass, so a
  usable GPU adapter exists in this environment** (upload/download round-trips,
  strided-elementwise-vs-CPU, sparse spmv/spmm all green). The source repo is not
  broken. **Remaining:** the git-dep *version-alignment* skew (hephaestus uses local
  path deps to the leto/mnemosyne/themis cluster) means Helios must consume it via a
  local `[patch]`/path (synchronized checkout), same wiring as G-11.
  *Update (H-010): CLOSED.* `helios-gpu` dispatches a real GPU kernel —
  `beam_transmission_into` computes `exp(-τ)` on the GPU (hephaestus-wgpu
  `NegOp`+`ExpOp`); a differential test vs CPU `f32::exp` passes on the live adapter.
  Wiring: replicated hephaestus's mnemosyne/moirai/hermes `[patch]` set so the
  leto→mnemosyne(git 1e014d25)→themis ^0.8 skew resolves to the local consistent
  cluster; hephaestus-wgpu consumed with default features (its `linalg` uses
  `leto-ops` ungated). *Evidence tier: verified — Helios GPU kernel runs + matches
  CPU (67 tests).* Remaining: throughput benchmark vs VoLO pending. (H-010b fused HU→μ
  kernel delivered 2026-07-02 via the hephaestus ADR-0004 authored-kernel
  seam — consumer-side `GpuAttenuationMapper`, differential vs solver 9/9.)

### Testing / tooling

- **G-8 (coverage):** No `cargo-llvm-cov` run yet; >80% core-logic coverage target
  unmeasured. Applies from first `[minor]`. → measure after H-003.

## Closed gaps

- **G-9 (CI):** CI now runs format, warning-denied Clippy, configured Nextest,
  doctests, warning-clean rustdoc, RustSec audit, and cargo-deny
  license/source policy. Benchmark classification is separately owned by the
  exact Atlas gate recorded under G-27.
- **G-6 (build hygiene):** Helios inherits the shared `D:/atlas/target` build dir
  via `repos/.cargo/config.toml`; no per-repo target. Verified this session.

## Integrated-workflow status

- **Integrated imaging-delivery workflow (clinical-realism gate): demonstrated end-to-end
  on synthetic/self-consistent data (H-041).** `helios-simulation/tests/end_to_end.rs`
  runs a shared μ through both branches (Radon→FBP→registration; helical MLC delivery→
  divergent-fan dose→scatter→DVH/gamma) with self-consistency oracles. What remains for
  the *clinical* therapy gate is orthogonal and environment-blocked: a licensed CT/plan
  dataset and an external MC reference engine (G-16), plus the anisotropic CC kernel
  (H-020g). The workflow *plumbing* across all layers is verified.

## Concurrent-agent status

- **RESOLVED (next cycle).** The peer `mnemosyne-arena` breakage that blocked `helios-gpu`
  and the moirai consumption last cycle has been reconciled by the peer (new mnemosyne
  commits landed; the crate compiles). The full workspace builds green again, and the
  moirai consumption (H-021b) was re-applied and verified — `simulate_helical_sinogram`
  dispatches per-projection work through moirai's `Adaptive` policy. Kept per the
  concurrent-agent discipline: the peer crate was never touched; the change was designed,
  reverted to stay green while blocked, then re-landed once upstream compiled.

## Residual risk register

- **G-20 (H-011b table interpolation).** The embedded dry-air, liquid-water,
  and cortical-bone `μ/ρ` values are transcribed from the NIST X-Ray Mass
  Attenuation Coefficients tables at their common 10 keV–20 MeV knots. The
  range excludes the selected tables' absorption-edge rows, and Helios defines
  between-knot behavior as native-precision log-linear interpolation. This is
  deliberately not represented as XCOM output: NIST documents log-log cubic
  spline fitting and explicit edge handling for XCOM. Exact table knots,
  boundaries, and the interpolation identity are value-semantic tests; an
  independent clinical-spectrum or Monte-Carlo validation remains outside this
  data-loading slice. *Evidence tier: source-value and analytical-contract
  tests.*

- Atlas upstream APIs may drift (multi-repo co-evolution); Helios pins the local
  synchronized checkout via `[patch]` and commits `Cargo.lock`. `ritk-dicom` is now
  consumed (H-004b) and is **skew-free** (no leto/mnemosyne/themis/eunomia cluster —
  only anyhow/arrayvec plus the RITK DICOM provider and `ritk-codecs`), so it needed no
  patch-cluster work. Helios has no direct dicom-rs dependency; the provider's
  parser implementation remains upstream-owned.
  Remaining ritk surfaces (`ritk-registration`) pull the burn stack and are heavier
  (G-5); add cross-repo contract tests as each is consumed.
- **G-18 — RESOLVED (H-043b).** The residency step landed: hephaestus gained a volume
  ray-integral kernel (`ray_line_integrals`) and `helios_gpu::GpuProjector` keeps μ
  on-device, projecting whole sinograms per dispatch. Measured 171×/371× vs the
  single-thread CPU projector (report `2026-07-02-gpu-projection-throughput.md`);
  differential per-ray agreement within a derived 1e-3 f32 bound. The *elementwise*
  `exp(−τ)` path remains transfer-bound by physics (documented; use the resident
  pipeline instead). "VoLO-competitive" is still unclaimable (no reference engine).
- **(historical) G-18 (performance, GPU transfer-bound).** The GPU-vs-CPU study (H-043,
  `validation_reports/2026-07-01-gpu-transmission-throughput.md`) shows the isolated
  `beam_transmission_into` kernel is memory-/transfer-bound: even on an RTX 5080 it
  reaches only ~0.5–0.72× a single-threaded CPU loop, because every call round-trips the
  buffer over PCIe for ~1 flop/element. This is a correct roofline result, not a defect.
  GPU throughput requires an **on-device fused pipeline** (H-043b) that keeps τ resident
  across HU→μ / projection / transmission so one CT upload + one sinogram download
  amortize many kernels. Until then the GPU path is a differentially-correct reference,
  not a speedup. The "competitive with VoLO-class throughput" gate additionally needs an
  external VoLO reference not available here. *Evidence tier: empirical (criterion, this
  machine).*
- **G-17 (tooling, coverage gate — link unblocked, attribution still empty).** Refined:
  the original blocker (the mingw `ld` bfd linker failing on `__llvm_profile_runtime` /
  profiler builtins) **is fixable** — `RUSTFLAGS="-Clink-arg=-fuse-ld=lld"` (LLVM `lld`,
  present in the MSYS2 ucrt64 toolchain) links the instrumented binaries, and the full
  suite runs under instrumentation (183 tests pass, 356 `.profraw` generated;
  `LLVM_COV`/`LLVM_PROFDATA` point at the MSYS2 llvm-cov 22.1.4 ≈ rustc-LLVM 22.1.3).
  A *distinct* secondary issue remains and is now **conclusively diagnosed**: source/
  region *attribution* is broken on this GNU target — `cargo llvm-cov report` gives
  **0 regions** and `grcov 0.10.5` (which uses its own profraw parser) gives an empty
  file table / **NaN%**, from the *same* 145 profraw. Two independent tools failing
  identically confirms the coverage-map is not read from the mingw (`x86_64-pc-windows-
  gnu`) instrumented binaries — a toolchain-level limitation, not a tool bug. Coverage %
  is therefore **not obtainable on this host**; it requires `x86_64-pc-windows-msvc` or a
  Linux CI container (H-060 re-scoped to CI). Test breadth is high (189 value-semantic
  tests across the CPU crates) but the coverage number is unquantified — not fabricated.
- Physical constants (G-2) are CODATA-2018/ICRU-90 values verified by inter-constant
  derivation tests, not by an external authoritative fetch this session; values are
  standard and cross-checked, but a future audit should confirm against the live
  NIST database.
- **G-16 (dose model fidelity, H-020d/H-020e).** *Partially addressed.* Stage 1
  (`deposit_ray_terma`/`accumulate_delivered_dose`) deposits primary terma along
  **parallel** beamlets; stage 2 (`scatter_superposition`, H-020e) now spreads it with
  a **separable-isotropic** deposition kernel, so lateral penumbra and depth build-up
  are present and energy-conserving (verified). The beam geometry now supports a
  divergent point-source fan (H-020f, `BeamGeometry::PointSource`; verified parallel
  limit + multi-row divergence). Still approximate vs a clinical collapsed-cone dose:
  inverse-square fluence falloff along the divergent fan is now modelled (H-020g,
  `deposit_ray_terma_diverging`; verified SAD→∞ limit + entry/exit steepening). The
  remaining approximation is the scatter kernel: separable-isotropic, not the
  anisotropic forward-peaked beam-aligned CC kernel, tracked as H-020h. Sufficient to exercise DVH/gamma
  on self-consistent phantoms; the therapy gamma/DVH clinical-agreement gate still
  needs the H-020g kernel upgrade AND a licensed real CT dataset AND an external
  Monte-Carlo/reference dose engine (VoLO/TOPAS/GATE/EGSnrc) — the last of which is
  **not runnable in this environment**, so that specific gate cannot be closed here and
  will not be fabricated. Evidence tier: analytical oracles (conservation, identity
  differential, symmetry); NOT validated against a reference dose engine.
