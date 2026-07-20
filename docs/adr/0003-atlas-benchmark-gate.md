# ADR 0003: Adopt the Atlas Criterion regression gate

- Status: Accepted
- Date: 2026-07-20
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

Helios pins Atlas merge `9bfb722` for both path-dependency checkout and
`tools/criterion-regression`. Pull-request CI:

1. checks out the pull-request base and candidate revisions on one runner;
2. copies the candidate benchmark sources into the baseline checkout so both
   revisions use one measurement instrument;
3. runs ABBA followed by its BAAB phase reversal;
4. retains the four Criterion comparison roots;
5. derives the per-case confidence from the complete benchmark family; and
6. delegates classification to Atlas
   `check-replicated-counterbalanced`.

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

The provider graph argument is independently pinned to Atlas `afd5e16`, whose
gitlinks match the Aequitas, Proteus, Hephaestus, Gaia, and Leto manifests
represented by `Cargo.lock`. The checkout action and Criterion implementation
remain pinned to their originating Atlas merge `9bfb722`; separating the graph
argument prevents tool provenance from becoming an obsolete provider snapshot.
Before measurement, CI resolves the historical baseline lock once against that
exact Ubuntu provider graph. Every measured baseline and candidate run then
uses `--locked`, and the delivered candidate lock is never regenerated.

## Rejected alternatives

- Keeping the same-run check remains tautological.
- A fixed percentage threshold discards the measured uncertainty.
- A Helios-owned Rust port would preserve duplicate statistical ownership.
- One ABBA block remains exposed to run-phase effects already falsified by the
  Apollo hosted canary recorded in Atlas ADR 0024.

## Consequences

- The benchmark classifier is an exact-revision CI dependency, not a Helios
  runtime dependency.
- Push events do not classify performance because they have no pull-request
  base contract.
- Candidate benchmark sources must compile against the baseline production
  revision. An incompatible instrument change fails visibly rather than
  producing a mixed-instrument claim.
- Static and synthetic evidence verifies classifier integration; only the
  hosted base/candidate lane supplies performance evidence.

## Verification

- Parse the workflow as YAML and scan it for the exact Atlas pin, four report
  roots, Nextest, doctests, and absence of the Python classifier.
- Run workspace format, warning-denied Clippy, configured Nextest, doctests,
  and warning-clean rustdoc locally.
- Require the exact pull-request head's Rust and benchmark jobs to pass.

## References

- [Atlas ADR 0024 at the pinned merge](https://github.com/ryancinsight/atlas/blob/9bfb722367a6c3911409d6b4619701c549b6d415/docs/adr/0024-criterion-regression-gate.md)
