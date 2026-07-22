# Chapter 16 — Gamma Index and Plan Verification

The gamma index is the standard plan QA metric for 3D dose comparison:

`	ext
γ(r_ref) = min_{r_eval} √[ (D_ref − D_eval)²/(δD)² + |r_ref − r_eval|²/(δr)² ]
`

A point passes if γ < 1.

`
ust
use helios_analysis::{gamma_index_3d, gamma_pass_rate};

let gamma = gamma_index_3d(
    &dose,        // evaluated distribution
    &reference,   // reference distribution
    0.03,         // ΔD tolerance (3 % of prescription)
    2.0,          // Δr distance tolerance (mm)
    5.0,          // dose threshold (ignore voxels < 5 Gy)
);

let pass_rate = gamma_pass_rate(&gamma, 1.0);
println!("3%/2mm pass rate: {:.1}%", pass_rate * 100.0);
`

## Self-Consistency Test

The tomotherapy workflow runs a **self-gamma** — the computed dose compared
against itself as reference. A perfect self-gamma always yields 100% pass
rate (within floating-point noise), confirming the metric implementation.

## Clinical Standard

AAPM TG-218 recommends 3%/2 mm criteria with 10% dose threshold.
Pass rate ≥ 95% is the typical clinical acceptance criterion.

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Dose-Volume Histograms](planning_dvh.md)