# HELIOS-ROI-VARIANCE-20261108 [patch]

Target: helios-analysis 0.1.1. Owner: Atlas scheduled correctness review.
Scope: shared intensity/dose ROI statistics; no provider, GPU, owner-branch,
registry, workflow, deployment or access changes.

## Specification and reproduction

For two samples [a, a+d], mean = a+d/2 and population standard deviation =
|d|/2. Adding a constant cannot change variance in exact arithmetic.
At base 4f90e77b8bd6baa5ae079fff0ce7a1c53347877c, native public-API tests
return zero standard deviation for both binary-exact witnesses:
- single precision [1024, 1024.25], expected 0.125;
- double precision [134217728, 134217730], expected 1.
The raw-moment subtraction E[x²] - E[x]² loses the small residual.

## Repair and verification

- [x] Replace raw moments with Welford's one-pass centered recurrence.
- [x] Keep one generic implementation and the caller's native precision.
- [x] Cover positive/negative/zero offsets, both public wrappers, and
      empty, reversed, singleton and partially out-of-bounds regions.
- [x] Native analysis suite: 40 passed, 0 skipped before version-only changes.
- [ ] Final revision checks, independent review and normal delivery hooks.
- [ ] Hosted gates, required approvals, provider qualification and deployment
      settings must pass before merge. Existing provider-cohort PR #183,
      hook PR #185 and protected rescue PR #173 remain untouched.

Residual scope: this closes cancellation from subtracting raw moments;
finite-range overflow and general correctly-rounded statistics are not proved.
No package publication or full Helios validation is claimed.
