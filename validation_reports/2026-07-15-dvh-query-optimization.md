# DVH threshold-query optimization — baseline comparison (H-062)

**Date:** 2026-07-15 · **Benchmark:** `helios-analysis/benches/dvh_queries.rs`
(Criterion, 100 samples) · **Machine:** Intel Core Ultra 9 285K, 24 logical
processors, single-threaded benchmark process.

## Workload

The benchmark builds one fixed 64³ `f64` volume, constructs one `Dvh`, and
issues 1,024 deterministic threshold queries. The scan arm counts qualifying
values directly from `Dvh::dose_sample()`; the production arm calls
`Dvh::volume_fraction_at_dose`. Both arms retain their accumulated result.

## Measured

| Arm | Median | 95% interval |
|---|---:|---:|
| Scan reference | 30.090 ms | 29.717–30.472 ms |
| Production binary bound | 29.229 μs | 28.426–30.023 μs |

The paired median ratio is **1,029×** for this repeated-query workload. A
pre-change scan run in the same checkout measured 28.820 ms [28.572, 29.076]
ms; the post-change scan control was 30.090 ms [29.717, 30.472] ms, so the
optimization claim uses the post-change control arm rather than attributing
machine variance to the production change.

## Invariant and residual risk

`Dvh::from_volume_masked` already stores doses sorted ascending with
`total_cmp`. For samples without NaN, the lower-bound predicate `value < dose`
is monotone and returns exactly the count satisfying `value >= dose`. A NaN
threshold returns zero and a sample containing NaN uses the previous direct
filter, preserving value semantics for non-finite input. The fallback is
intentionally O(n) because unordered values cannot satisfy the binary-search
monotonicity contract; valid finite dose plans take the O(log n) path.

No allocation or ownership transfer was added to the query path. Focused
nextest passes 34/34, including finite boundary and NaN differential cases.
