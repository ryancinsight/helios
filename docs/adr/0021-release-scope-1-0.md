# ADR 0021: Helios release scope for 1.0

Status: Accepted  \
Date: 2026-11-15  \
Change class: [arch]  \
Board item: HELIOS-RELEASE-LADDER

## Context

Helios is the stack's radiative-transport and dose owner (photon transport,
CT imaging, dose accumulation, inverse planning). The strategic board carries
its history as done slices plus fifteen open items: dose-engine validation
and planning wiring, projection specification, DICOM-RT interop, the Python
surface, book re-grounding, and board/lock hygiene. No milestone ladder
orders them toward a stable 1.0, and no record states what 1.0 freezes.

## Decision

Adopt a six-rung ladder to 1.0. The ladder opens `backlog.md`; every open
item carries exactly one `milestone:` tag.

- **M-dose-engine**: delivered dose carries evidence above analytical
  self-consistency (H-112), and inverse planning drives the dose engine
  through a `DoseInfluence` producer (H-114). Exit: the reference-engine or
  published-benchmark comparison is green within derived tolerances, or the
  absence is a recorded bounded decision; the producer wires planning to
  simulation.
- **M-imaging**: the projector is exact and provider-delegated (H-011d,
  H-021c, H-113). Exit: ADR 0019 cases 1-8 pass in f32/f64, fan rays match
  the parallel projector within the roundoff bound, and no declared Atlas
  provider stays locally substituted.
- **M-interop**: clinical and scripting boundaries are typed (H-115, H-118,
  H-118a) with the charset dependency externally blocked (H-073, reopens on
  a maintained `dicom-encoding` release). Exit: DICOM-RT I/O round-trips
  structure sets, the Python surface is typed with PEP 561 metadata, and
  H-073 stays a recorded external hold.
- **M-book**: the book's claims resolve against the tree (H-110, H-111,
  H-116). Exit: the scripted chapter-body check reports zero unresolved
  paths, `mdbook test` executes a non-zero test count, and the figure tree
  matches the SUMMARY.
- **M-board**: hygiene and locks are green (H-117, H-120, H-121). Exit: no
  stale statuses or duplicate IDs, the benchmark runner is partitioned
  within the file target, and the locked workspace compiles against one
  coherent provider graph.
- **M-release**: stable transport-and-dose release (**authority-gated**). Exit:
  `cargo semver-checks`, publish dry-run, `cargo doc`, the full
  nextest/doctest suite green, and zero open correctness items.
  **Veto point:** release authority sets the API-freeze contents; this
  milestone never forces `publish = true` without it.
- **M-ecosystem**: consumer-owned coupling (harmonia therapy-response
  coupling, asclepius bioheat consumers). Exit: consumer receipts on their
  own boards; does not gate 1.0.

### Veto points

1. The 1.0 API-freeze contents (release authority sets them).
2. The H-112 evidence tier: a reference-engine comparison or a recorded,
   bounded decision with the exact blocker; this ADR does not pre-pick it.

## Rejected alternatives

- Freezing 1.0 on analytical self-consistency alone: H-112 records that
  dose without an independent tier repeats the oracle gap the stack's other
  members closed with published benchmarks.
- Counting H-073 as a 1.0 gate: the charset dependency is external with no
  maintained replacement published; it stays a recorded hold outside the
  ladder's blocking path.

## Overturning evidence

A consumer failure traceable to a frozen API, or a dose-validation
mismatch against the chosen reference, reopens the affected rung.

## Revision

2026-11-15: the ladder block `M-1.0` is `M-release`. The board lint item
pattern accepts a hyphen segment of `[A-Z0-9]+`, so it read `M-1` as an item
ID because the segment `1` is a digit; `M-release` has no such segment.
