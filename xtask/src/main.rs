use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::{Path, PathBuf};

mod check_figures;
mod migration_audit;
mod prebook;

#[derive(Parser)]
#[command(name = "xtask")]
#[command(about = "Helios automation tasks")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Audit legacy nalgebra/ndarray/burn/tokio/rayon migration surface.
    LegacyMigrationAudit,
    /// Refresh the legacy migration allowlist baseline file.
    RefreshLegacyAllowlist,
    /// Audit direct Burn migration surface during the Coeus cleanup.
    BurnMigrationAudit,
    /// Refresh the Burn migration allowlist baseline file.
    RefreshBurnAllowlist,
    /// Verify the deterministic figure set committed at `docs/book/figures/`
    /// and write a `MANIFEST.json` next to it. The manifest hash is the
    /// byte fingerprint each file currently has on disk; re-running on
    /// unchanged inputs produces byte-identical output (CI evidence chain).
    Prebook,
    /// Verify SUMMARY.md + README.md figure links are listed in FIGURE_SPECS.
    /// Returns exit code 1 on drift so the CI gate fails loudly when the
    /// SSOT contract between figures and book source breaks.
    CheckFigures,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let root = workspace_root();

    match cli.command {
        Command::LegacyMigrationAudit => {
            migration_audit::print_legacy_migration_audit(&root)
        }
        Command::RefreshLegacyAllowlist => {
            migration_audit::refresh_legacy_allowlist(&root)
        }
        Command::BurnMigrationAudit => {
            migration_audit::print_burn_migration_audit(&root)
        }
        Command::RefreshBurnAllowlist => migration_audit::refresh_burn_allowlist(&root),
        Command::Prebook => run_prebook(&root),
        Command::CheckFigures => run_check_figures(&root),
    }
}

fn run_check_figures(workspace_root: &Path) -> Result<()> {
    let report = check_figures::run_check_figures(workspace_root)?;
    let exit_code = check_figures::print_check_figures_report(&report);
    if exit_code != 0 {
        std::process::exit(exit_code);
    }
    Ok(())
}

fn run_prebook(workspace_root: &Path) -> Result<()> {
    let report = prebook::run_prebook(workspace_root)?;
    println!(
        "prebook: verified {} figures under {}",
        report.entries.len(),
        workspace_root.join("docs/book/figures").display()
    );
    for entry in &report.entries {
        let source = match entry.source_example {
            Some(e) => format!("{}::{}", e.crate_name, e.example_name),
            None => "schematic".to_string(),
        };
        println!(
            "  - {:<32} {:>7}B  sha256:{}  ({})",
            entry.name, entry.bytes, entry.sha256_16, source
        );
    }
    println!("manifest: {}", report.manifest_path.display());
    Ok(())
}

fn workspace_root() -> PathBuf {
    let xtask_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    xtask_dir
        .parent()
        .map(PathBuf::from)
        .unwrap_or_else(|| xtask_dir)
}
