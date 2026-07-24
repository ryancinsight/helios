//! Prebook / figure-management for the Helios mdbook.
//!
//! Build the deterministic figure set committed at `docs/book/figures/`
//! and emit a `MANIFEST.json` manifest that downstream tooling (link-checker
//! pre-flight, CI evidence chain) can fingerprint.
//!
//! # Contract
//!
//! - `FIGURE_SPECS` is the **single source of truth** (SSOT) for the
//!   figure set. `helios/docs/book/SUMMARY.md` and `helios/docs/book/README.md`
//!   stay manually consistent with this list.
//! - Each figure is `HandAuthored`: the `.svg` is committed by hand to
//!   mirror the deterministic output of a specific example.  Both the SVG
//!   and the example's printed values are deterministic, so two `prebook`
//!   invocations on unchanged inputs produce byte-identical output.
//! - The manifest hash is the first 16 hex chars of SHA-256 over the
//!   figure file content; 64 bits is plenty for a book-scale figure set
//!   (collision probability for 7 entries ≈ 2.7 × 10⁻¹⁵).

use anyhow::{anyhow, Context, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

/// A figure referenced by the Helios mdbook.
///
/// `source_example` is `Some(_)` for figures that mirror a runnable
/// example's deterministic stdout / staged PNG; `None` for purely
/// authored schematics (workflow / architecture diagrams).
#[derive(Debug, Clone, Copy, Serialize)]
pub struct FigureSpec {
    /// File name (relative to `docs/book/figures/`).
    pub name: &'static str,
    /// Optional runnable example this figure mirrors.
    pub source_example: Option<ExampleRef>,
    /// One-line summary used in the README figure index.
    pub summary: &'static str,
}

/// A runnable example reference: `cargo run -p <crate_name> --example <example_name>`.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct ExampleRef {
    pub crate_name: &'static str,
    pub example_name: &'static str,
}

/// Authoritative figure list — kept synchronised with `SUMMARY.md` and
/// `README.md` figure indexes by hand.
pub const FIGURE_SPECS: &[FigureSpec] = &[
    FigureSpec {
        name: "photon_attenuation_depth.svg",
        source_example: Some(ExampleRef {
            crate_name: "helios-physics",
            example_name: "photon_attenuation",
        }),
        summary: "Exponential transmission T(x) vs depth in water at 100 keV (mu = 0.171 cm^-1, depths 0..20 cm).",
    },
    FigureSpec {
        name: "ct_calibration_curve.svg",
        source_example: Some(ExampleRef {
            crate_name: "helios-physics",
            example_name: "photon_attenuation",
        }),
        summary: "Hounsfield-unit vs relative-electron-density calibration across six reference materials.",
    },
    FigureSpec {
        name: "radon_sinogram_disk.svg",
        source_example: Some(ExampleRef {
            crate_name: "helios-imaging",
            example_name: "radon_sinogram",
        }),
        summary: "Single-angle (theta = 0) chord profile for a 30 mm radius disk phantom.",
    },
    FigureSpec {
        name: "dvh_curve.svg",
        source_example: Some(ExampleRef {
            crate_name: "helios-analysis",
            example_name: "dvh_analysis",
        }),
        summary: "Differential DVH of the synthetic Gaussian dose phantom used in dvh_analysis.",
    },
    FigureSpec {
        name: "dose_slice_heatmap.svg",
        source_example: Some(ExampleRef {
            crate_name: "helios-simulation",
            example_name: "tomotherapy_workflow",
        }),
        summary: "Central-slice dose heatmap of the helical TomoTherapy workflow output.",
    },
    FigureSpec {
        name: "helical_mlc_fluence.svg",
        source_example: Some(ExampleRef {
            crate_name: "helios-simulation",
            example_name: "tomotherapy_workflow",
        }),
        summary: "21-leaf x 40-projection MLC leaf-open-time sinogram for the target aperture.",
    },
    FigureSpec {
        name: "architecture_stack.svg",
        source_example: None,
        summary: "Six-layer top-down Helios -> Atlas stack architecture diagram.",
    },
];

/// Single manifest entry — produced by [`run_prebook`] and serialised as
/// JSON to `docs/book/figures/MANIFEST.json`. Byte-deterministic across
/// runs of identical inputs (no timestamps, no ordering dependence).
#[derive(Debug, Clone, Serialize)]
pub struct ManifestEntry {
    pub name: String,
    pub source_example: Option<ExampleRef>,
    pub sha256_16: String,
    pub bytes: usize,
    pub summary: String,
}

/// Aggregated report returned to the CLI surface.
#[derive(Debug, Clone)]
pub struct PrebookReport {
    pub entries: Vec<ManifestEntry>,
    pub manifest_path: PathBuf,
}

/// Run prebook against `workspace_root`. For each spec, verifies the
/// figure file exists at `<workspace_root>/docs/book/figures/<name>`,
/// hashes it, and writes a deterministic `MANIFEST.json` next to the
/// figures. The returned report lists every verified entry.
pub fn run_prebook(workspace_root: &Path) -> Result<PrebookReport> {
    let figs_dir = workspace_root.join("docs/book/figures");
    if !figs_dir.is_dir() {
        return Err(anyhow!(
            "figures directory not found: {} (prebook expects <workspace>/docs/book/figures/)",
            figs_dir.display()
        ));
    }
    validate_figure_files(&figs_dir)?;

    // Iterate specs in declaration order — stable, JSON-ordered manifest.
    let mut entries: Vec<ManifestEntry> = Vec::with_capacity(FIGURE_SPECS.len());
    for spec in FIGURE_SPECS {
        let path = figs_dir.join(spec.name);
        let bytes =
            fs::read(&path).with_context(|| format!("reading figure file {}", path.display()))?;
        let sha = sha256_hex_first_16(&bytes);
        entries.push(ManifestEntry {
            name: spec.name.to_owned(),
            source_example: spec.source_example,
            sha256_16: sha,
            bytes: bytes.len(),
            summary: spec.summary.to_owned(),
        });
    }

    // Serialise with sorted keys (serde_json default) + no pretty-print
    // trailing whitespace; deterministic across runs and machines.
    let manifest_path = figs_dir.join("MANIFEST.json");
    let json = serde_json::to_string(&entries).context("serialising MANIFEST.json entries")?;
    fs::write(&manifest_path, format!("{json}\n"))
        .with_context(|| format!("writing manifest {}", manifest_path.display()))?;

    Ok(PrebookReport {
        entries,
        manifest_path,
    })
}

/// Verify that the committed SVG set and [`FIGURE_SPECS`] are identical.
///
/// Missing files are reported by [`run_prebook`] while hashing the declared
/// entries; this check closes the opposite direction so an unlisted SVG
/// cannot silently bypass the manifest SSOT.
pub fn validate_figure_files(figs_dir: &Path) -> Result<()> {
    let declared: BTreeSet<&str> = FIGURE_SPECS.iter().map(|spec| spec.name).collect();
    let actual: BTreeSet<String> = fs::read_dir(figs_dir)
        .with_context(|| format!("reading figure directory {}", figs_dir.display()))?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .filter(|path| path.extension().is_some_and(|ext| ext == "svg"))
        .filter_map(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .map(str::to_owned)
        })
        .collect();

    let extras: Vec<&str> = actual
        .iter()
        .filter(|name| !declared.contains(name.as_str()))
        .map(String::as_str)
        .collect();
    if !extras.is_empty() {
        return Err(anyhow!(
            "unlisted SVG figure(s) in {}: {}",
            figs_dir.display(),
            extras.join(", ")
        ));
    }

    Ok(())
}

/// SHA-256 hex digest; only the first 16 hex chars are surfaced in the
/// manifest (24 hex chars would still collide in ~10^19 entries; 16
/// hex chars at ~6.5·10^7 distinct values is plenty for a book-scale
/// figure set and the file content is small + committed, so anyone
/// needing the full hash can recompute).
pub fn sha256_hex_first_16(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let digest = hasher.finalize();
    let hex = format!("{:x}", digest);
    hex[..16].to_owned()
}
