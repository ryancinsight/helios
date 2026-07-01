# Helios (Helios-rs)

A modular Rust workspace for **unified radiation-therapy simulation/planning and
radiation imaging**, built natively on the [Atlas](../..) technology stack.

Helios targets VoLO-class GPU-accelerated TomoTherapy — helical delivery with
binary-MLC modulation, couch translation, and MVCT guidance — alongside
high-fidelity radiation-imaging simulation (MVCT forward projection,
reconstruction, on-board kV/MV imaging, and IGRT workflows).

> **Status:** Sprint 1 (Foundation). `helios-core` is implemented and verified; the
> higher layers are being built out per the sprint roadmap. See
> [`CHECKLIST.md`](CHECKLIST.md), [`backlog.md`](backlog.md), and
> [`gap_audit.md`](gap_audit.md).

## Architecture at a glance

```
helios-python                         PyO3 bindings (thin)
   │
helios-planning  helios-simulation  helios-imaging     application
   │                │                 │
helios-analysis   helios-solver ── helios-gpu           compute
   │                │
helios-physics ─────┘
   │
helios-domain                                           domain
   │
helios-math                                             numerics (Scalar seam)
   │
helios-core                                             foundation
```

Strictly unidirectional layering; Atlas crates (ritk, gaia, hephaestus, moirai,
coeus, consus, leto, hermes, mnemosyne, themis, apollo) are consumed as remote git
dependencies. Full crate responsibilities and the Atlas dependency map are in
[`ARCHITECTURE.md`](ARCHITECTURE.md).

## Building

```sh
# Rust stable (edition 2021); toolchain pinned via rust-toolchain.toml.
cargo build
cargo clippy --all-targets --all-features -- -D warnings
cargo nextest run          # test-time budget: 30 s slow / 60 s terminate
cargo test --doc           # doctests (nextest does not run these)
```

The workspace uses the shared Atlas test-time budget: a test that crosses 30 s is a
performance defect to profile and optimize, never a limit to raise.

## Sprint roadmap

1. **Sprint 1 — Foundation:** workspace skeleton, `helios-core`, `helios-math`
   seam, `helios-domain` + ritk/gaia integration (CT/MVCT + geometry). *(in
   progress)*
2. **Sprint 2 — GPU foundation:** hephaestus+moirai; first deterministic dose +
   imaging projection kernels in `helios-solver`/`helios-gpu`.
3. **Sprint 3 — Delivery:** binary-MLC modeling + full helical TomoTherapy delivery
   kinematics in `helios-simulation` with basic MVCT.
4. **Sprint 4 — Planning & imaging:** `helios-planning` (coeus) + dedicated
   `helios-imaging`.
5. **Sprint 5 — End-to-end:** workflow validation, Python bindings, performance
   benchmarking.

## Validation targets

- **Therapy:** gamma analysis (3%/2 mm), DVH agreement vs VoLO, and reference Monte
  Carlo (TOPAS, GATE, EGSnrc) on phantom/clinical cases.
- **Imaging:** MVCT reconstruction accuracy, noise, contrast, spatial resolution vs
  published TomoTherapy MVCT data.
- **Performance:** GPU scaling and timing competitive with VoLO-class throughput.
- **Software:** zero Clippy warnings on production paths, >80% core coverage,
  property-based testing, benchmarks with recorded baselines.

## License

Dual-licensed under MIT or Apache-2.0.
