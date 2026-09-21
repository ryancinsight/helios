# ADR 0003: Benchmarks are a local instrument; CI keeps only the bench smoke

- Status: Accepted
- Date: 2026-07-20
- Revision: 2026-09-21 — the hosted timing gate is removed. The
  single-runner `benchmark-regression` job owes eight measured legs
  (four `compare_pair` calls over two revisions) against a 60-minute
  cap and dies red mid-schedule (exit 124, recorded on the
  `ci/bench-pair-budget` branch against PR #88's schedule); splitting
  the same legs across four hosted runners (draft PR #90) keeps the
  defect and adds cross-runner variance. Shared-runner wall-clock
  timings are noise, not evidence, and the Atlas-owned classifier
  contract (`tools/criterion-regression/README.md`) already states
  that the four pairs run within the committed local timing budget
  while CI only smoke-runs benchmarks. The decision below records the
  correction. What follows is the now-valid contract; git history is
  the archive of the superseded hosted-gate text.
- Class: `[arch]` `[patch]`

## Context

Helios CI ran the candidate benchmark suite once, serialized that Criterion
tree through a copied Python script, and immediately compared the same tree
with itself. The gate could not detect a regression, and its empirical 15%
threshold had no error model. The Rust job also bypassed the committed Nextest
timeout profile with bare `cargo test`.

Atlas ADR 0024 owns the statistical and cross-repository contract. Helios needs
only the consumer orchestration for its four Criterion benchmark targets.

## Decision

Helios pins Atlas merge `9bfb722` for `tools/criterion-regression`.
Benchmark timing evidence is produced only by the committed local
instrument; pull-request CI never classifies performance:

1. the operator runs the paired schedule locally via the committed
   runner (`xtask bench-replicated`: `A B B A` followed by its `B A A B`
   phase reversal) on one controlled host, holding the candidate
   benchmark sources constant across both revisions;
2. each leg runs the four declared Criterion binaries
   (`helios-analysis:dvh_queries`, `helios-gpu:projection_throughput`,
   `helios-gpu:transmission_throughput`,
   `helios-solver:scatter_superposition`) with abbreviated sampling sized
   to the derived 1500 s suite bound (measured 1101 s same-revision
   calibration on the reference host, floored by ~15 s single iterations
   in `projection_throughput`);
3. the four retained Criterion comparison roots plus the derived
   confidence are attached to the PR as the performance evidence, and
   classification delegates to Atlas `check-replicated-counterbalanced`;
4. CI runs only the single-iteration bench smoke
   (`cargo test --benches`, equivalently Criterion `--test`) inside the
   standard test budget; a smoke failure blocks merge, a timing
   comparison never does;
5. the `benchmark-regression` and `classify` CI jobs are deleted, and
   draft PR #90 (parallel hosted pairs) closes with this ADR as its
   verdict.

The complete schedule is `A B B A B A A B`. Baseline and candidate each occupy
positions with sum 18 and squared sum 102, balancing exposure to constant,
linear, and quadratic period terms. Atlas requires a slowdown to reproduce
inside and across both blocks, controls family-wise false regressions at 5%,
and fails closed on incomplete or mismatched evidence.

The Rust job installs pinned cargo-nextest, cargo-audit, and cargo-deny
versions; runs the committed `ci` profile; runs doctests separately; and
enforces RustSec, license, and dependency-source policy. Only
RUSTSEC-2021-0153 is quarantined: current `dicom-encoding` requires the
unmaintained charset crate unconditionally, and the advisory reports no known
vulnerability. The copied Python classifier is deleted.

The binding job builds the abi3 extension with pinned Maturin and executes the
value-semantic Pytest suite against the installed wheel. This keeps Python as a
tested FFI boundary over the Rust cores rather than an unverified packaging
artifact.

The provider graph and checkout action are pinned together at Atlas merge
`05b7f5d`. Its gitlinks include public Asclepius `ceb8b6d` and the Hephaestus
revision that shares Helios's Aequitas identity, matching the Proteus, Gaia,
and Leto manifests represented by `Cargo.lock`. The Criterion implementation
remains pinned to its audited Atlas merge `9bfb722`; provider advancement does
not change the measurement instrument. Before measurement, CI resolves the
historical baseline lock once against that exact Ubuntu provider graph. Every
measured baseline and candidate run then uses `--locked`, and the delivered
candidate lock is never regenerated.

## Rejected alternatives

- Keeping the same-run check remains tautological.
- A fixed percentage threshold discards the measured uncertainty.
- A Helios-owned Rust port would preserve duplicate statistical ownership.
- One ABBA block remains exposed to run-phase effects already falsified by the
  Apollo hosted canary recorded in Atlas ADR 0024.
- Splitting the eight hosted legs across four runners (draft PR #90) keeps
  the category error and adds cross-runner variance: each leg still times
  on a shared runner, and the classifier's replication contract assumes
  co-located pairs on one controlled host.

## Consequences

- The benchmark classifier is an exact-revision CI dependency, not a Helios
  runtime dependency.
- Push events do not classify performance because they have no pull-request
  base contract.
- Candidate benchmark sources must compile against the baseline production
  revision. An incompatible instrument change fails visibly rather than
  producing a mixed-instrument claim.
- Static and synthetic evidence verifies classifier integration; only a
  controlled local host supplies performance evidence, never a shared CI
  runner.
- Each measured revision runs only the declared `harness = false` benchmark
  binaries. Workspace library targets remain part of the Rust correctness job
  and do not receive Criterion command-line arguments.
- A dependency-only candidate can change linked-code layout and therefore
  throughput even when the measured Rust source is unchanged. A replicated
  regression remains a production defect; the response is to optimize the hot
  path, never to weaken the instrument or statistical classifier.

## Verification

- Parse the workflow as YAML and confirm the `benchmark-regression` and
  `classify` jobs are absent while a bench-smoke step runs
  `cargo test --benches` inside the standard test budget.
- Run the committed local runner (`cargo xtask bench-replicated --help`
  documents the paired schedule); one full local schedule completes
  inside the derived 1500 s suite bound on the reference host class,
  with per-pair elapsed times printed for attribution.
- Run workspace format, warning-denied Clippy, configured Nextest, doctests,
  and warning-clean rustdoc locally.
- Pin the scatter-convolution rewrite against the unchanged Criterion
  instrument. Local paired evidence on the development host reports
  50.46% lower median time at 32³ and 51.02% at 64³; a bitwise differential
  test covers every axis and asymmetric boundary truncation.
- Require the exact pull-request head's Rust and bench-smoke jobs to pass;
  timing evidence arrives as PR-attached local reports, never as a CI gate.

## References

- [Atlas ADR 0024 at the pinned merge](https://github.com/ryancinsight/atlas/blob/9bfb722367a6c3911409d6b4619701c549b6d415/docs/adr/0024-criterion-regression-gate.md)
