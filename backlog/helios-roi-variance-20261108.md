<a id="helios-roi-variance-20261108"></a>
## HELIOS-ROI-VARIANCE-20261108 — Preserve translated ROI noise [patch] — blocked
- priority: correctness
- owner: Atlas scheduled correctness review; target: helios-analysis 0.1.1.
- scope: shared intensity/dose ROI statistics; provider and owner work untouched.
- Defect: raw moments erase nonzero variance. A drifting mean cancels adjacent-value variance; pre-normalization division erases subnormal center terms.
- Specification: a balanced sample of a and a+d has population standard deviation |d|/2.
- Repair: normalize before compensated center accumulation and a nonnegative residual norm; preserve the original returned mean and native precision.
- [x] Cover one-ULP, signed/reversed offsets, repeated 3-D samples, both APIs and empty/singleton/clipped regions.
- [x] Analysis/planning/simulation: 100 native passes, 0 skips; formatting, strict Clippy and Rustdoc pass; 0 doctests.
- [ ] Qualify the frozen revision independently and in hosted gates; PR creation stays stopped after HTTP 403.
- [ ] Qualify provider migration and deployment settings before merge; PRs #183, #185 and #173 remain untouched.
- Residual scope: the returned mean retains its historical range/rounding limits; no universal correctly-rounded statistics claim.
- No merge, publication or deployment is authorized by this checklist's review status.
