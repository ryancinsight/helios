use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

mod migration_audit;

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
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::LegacyMigrationAudit => {
            migration_audit::print_legacy_migration_audit(&workspace_root())
        }
        Command::RefreshLegacyAllowlist => {
            migration_audit::refresh_legacy_allowlist(&workspace_root())
        }
        Command::BurnMigrationAudit => {
            migration_audit::print_burn_migration_audit(&workspace_root())
        }
        Command::RefreshBurnAllowlist => migration_audit::refresh_burn_allowlist(&workspace_root()),
    }
}

fn workspace_root() -> PathBuf {
    let xtask_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    xtask_dir
        .parent()
        .map(PathBuf::from)
        .unwrap_or_else(|| xtask_dir)
}
