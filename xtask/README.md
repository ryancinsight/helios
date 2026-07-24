# xtask — helios mdbook automation

Build + validation automation for the `helios/docs/book/` mdbook. Sources:
`src/{prebook,check_figures,migration_audit}.rs`. Run via
`cargo run -p xtask -- <subcommand>` or `cargo xtask <subcommand>` (cargo
alias).

## `cargo xtask prebook`

Verify + SHA-256-hash the deterministic figure set committed at
`docs/book/figures/` into `MANIFEST.json`. The manifest is
byte-deterministic across repeated runs on unchanged inputs (CI evidence
chain). Source: `src/prebook.rs`.

## `cargo xtask check-figures`

SSOT drift verification between `docs/book/SUMMARY.md` +
`docs/book/README.md` figure links and the canonical `FIGURE_SPECS` list
in `src/prebook.rs`. Source: `src/check_figures.rs`.

### What it checks

For every markdown link of the form `[Title](figures/<name>.svg)` in
`SUMMARY.md` and `README.md`: the `<name>.svg` must appear in
`FIGURE_SPECS`. Conversely, every `FigureSpec` in `FIGURE_SPECS` must
have a matching `figures/<name>.svg` link in `SUMMARY.md` or `README.md`.

Drift in either direction is reported.

The parser is a pure byte-scan over `figures/` substrings (no regex dep
— keeps the `xtask` workspace graph lean). It terminates a figure name
on `)` or ASCII whitespace and accepts the token only if it ends with
`.svg` and contains no further `/`. Backtick code spans are not
specifically recognised, so put figure references in real markdown
links, not in inline code.

### When to run

- Before any PR that adds or removes a figure link in `SUMMARY.md` /
  `README.md`, or a `FigureSpec` entry in `src/prebook.rs`.
- As part of CI on every PR (recommended wiring; not yet committed).

### Exit codes

| Code | Meaning |
|------|---------|
| `0`  | `SSOT_IN_SYNC` — every docs figure link is listed in `FIGURE_SPECS`, and vice versa. |
| `1`  | Drift detected. Sub-reports: `DRIFT_DOCS_NOT_IN_SPECS` (a docs link has no matching spec) or `DRIFT_SPECS_NOT_IN_DOCS` (a spec has no matching docs link). |

### Output

```
check-figures: SSOT drift verification (SUMMARY.md + README.md vs FIGURE_SPECS)

  FIGURE_SPECS entries  : 7
  Docs references       : 7
  In sync (intersection): 7

SSOT_IN_SYNC: every docs figure link is listed in FIGURE_SPECS.
```

On drift, each orphan line is printed with the source file, line
number, and figure filename, e.g.:
`- SUMMARY.md:L42  (orphan_figure.svg)`.

## `FIGURE_SPECS` SSOT contract

`src/prebook.rs` declares `pub const FIGURE_SPECS: &[FigureSpec]`. Each
entry has:

```rust
pub struct FigureSpec {
    pub name: &'static str,                // e.g. "photon_attenuation_depth.svg"
    pub source_example: Option<ExampleRef>, // (crate_name, example_name) that generates the figure
    pub summary: &'static str,             // one-line description
}
```

`FIGURE_SPECS` is the single source of truth for the committed figure
set: every committed SVG under `docs/book/figures/` must be listed
here, and every entry here must be linked from `SUMMARY.md` or
`README.md`. The `check-figures` subcommand enforces that contract.

### Adding a new figure

1. Add a `FigureSpec` entry in `src/prebook.rs` `FIGURE_SPECS`.
2. Commit the SVG under `docs/book/figures/<name>.svg`.
3. Link it from `docs/book/SUMMARY.md` (chapter entry) or
   `docs/book/README.md` (figure index).
4. Run `cargo xtask check-figures` to confirm SSOT; run `cargo xtask
   prebook` to regenerate `MANIFEST.json` with the new SHA-256 entry.

## Other commands

- `cargo xtask migration-audit` — legacy surface audit (see
  `src/migration_audit.rs`).
- `cargo xtask refresh-legacy-allowlist` — regenerate
  `xtask/legacy_surface.allowlist`.

## Directory layout

```
helios/xtask/
├── Cargo.toml
├── README.md                       # this file
├── src/
│   ├── main.rs                     # subcommand enum + dispatch
│   ├── prebook.rs                  # FIGURE_SPECS + prebook + MANIFEST.json
│   ├── check_figures.rs            # SSOT drift lint
│   └── migration_audit.rs          # legacy surface scanner
├── burn_surface.allowlist          # tracked
└── legacy_surface.allowlist        # tracked
```
