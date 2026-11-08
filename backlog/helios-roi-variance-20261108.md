# HELIOS-ROI-VARIANCE-20261108 [patch]
Target: helios-analysis 0.1.1. Owner: Atlas scheduled correctness review.
Scope: shared intensity/dose ROI statistics; provider and owner work untouched.
Base: 4f90e77b8bd6baa5ae079fff0ce7a1c53347877c.
Defect: raw-moment cancellation reports zero noise for [1024, 1024.25] in single precision (expected 0.125) and [134217728, 134217730] in double precision (expected 1).
Specification: a balanced sample of a and a+d has population standard deviation |d|/2.
- [x] Corrected two-pass centered variance; original mean and native precision retained.
- [x] Cover one-ULP, signed/reversed offsets, repeated 3-D samples, both public APIs and empty/singleton/clipped regions.
- [x] Analysis/planning/simulation: 94 native tests pass, 0 skipped.
- [x] Formatting, strict Clippy and Rustdoc pass; doctest collection has 0 tests.
- [ ] Final independent review, normal hooks, hosted gates and required approvals.
- [ ] Qualify provider migration and deployment settings before merge; PRs #183, #185 and #173 remain untouched.
Residual scope: general finite-range overflow and correctly-rounded statistics are not proved.
No package publication or full Helios qualification is claimed.
