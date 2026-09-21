//! Revision checkout materialization and benchmark-source pinning.

use super::BenchReplicatedArgs;
use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Bench source trees held constant across both revisions.
const BENCH_DIRS: &[(&str, &str)] = &[
    ("helios-analysis", "crates/helios-analysis/benches"),
    ("helios-gpu", "crates/helios-gpu/benches"),
    ("helios-solver", "crates/helios-solver/benches"),
];

fn copy_dir(src: &Path, dst: &Path) -> Result<()> {
    for entry in std::fs::read_dir(src).with_context(|| format!("read {}", src.display()))? {
        let entry = entry?;
        let dst_path = dst.join(entry.file_name());
        if entry.file_type()?.is_dir() {
            std::fs::create_dir_all(&dst_path)?;
            copy_dir(&entry.path(), &dst_path)?;
        } else {
            std::fs::copy(entry.path(), &dst_path)?;
        }
    }
    Ok(())
}

struct TempRoot {
    path: PathBuf,
}

impl Drop for TempRoot {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

pub(super) struct Materialized {
    pub(super) baseline_dir: PathBuf,
    pub(super) candidate_dir: PathBuf,
    /// Owned temp root (ephemeral worktrees); `None` for explicit dirs,
    /// whose lifetime the operator owns.
    _temp_root: Option<TempRoot>,
    pub(super) run_root: PathBuf,
}

fn ensure_worktree(repo: &Path, rev: &str, path: &Path) -> Result<()> {
    if path.exists() {
        bail!(
            "worktree path {} already exists; remove it first",
            path.display()
        );
    }
    let status = Command::new("git")
        .arg("-C")
        .arg(repo)
        .arg("worktree")
        .arg("add")
        .arg("--detach")
        .arg(path)
        .arg(rev)
        .status()
        .context("git worktree add")?;
    if !status.success() {
        bail!("git worktree add {rev} failed with {status}");
    }
    Ok(())
}

fn remove_worktree(repo: &Path, path: &Path) {
    let _ = Command::new("git")
        .arg("-C")
        .arg(repo)
        .arg("worktree")
        .arg("remove")
        .arg("--force")
        .arg(path)
        .status();
}

pub(super) fn materialize(
    args: &BenchReplicatedArgs,
    workspace_root: &Path,
) -> Result<Materialized> {
    let run_root = workspace_root.join("target").join("bench-replicated");
    std::fs::create_dir_all(run_root.join("reports"))?;
    match (&args.baseline_dir, &args.candidate_dir) {
        (Some(base), Some(cand)) => Ok(Materialized {
            baseline_dir: base.clone(),
            candidate_dir: cand.clone(),
            _temp_root: None,
            run_root,
        }),
        (None, None) => {
            let root = std::env::temp_dir()
                .join(format!("helios-bench-replicated-{}", std::process::id()));
            if root.exists() {
                std::fs::remove_dir_all(&root)?;
            }
            std::fs::create_dir_all(&root)?;
            let baseline_dir = root.join("helios-baseline");
            let candidate_dir = root.join("helios-candidate");
            ensure_worktree(workspace_root, &args.baseline, &baseline_dir)?;
            if let Err(e) = ensure_worktree(workspace_root, &args.candidate, &candidate_dir) {
                remove_worktree(workspace_root, &baseline_dir);
                return Err(e);
            }
            Ok(Materialized {
                baseline_dir,
                candidate_dir,
                _temp_root: Some(TempRoot { path: root }),
                run_root,
            })
        }
        _ => bail!("--baseline-dir and --candidate-dir must be given together"),
    }
}

/// Hold the candidate benchmark instrument constant: the baseline measures
/// with the candidate's bench sources, so only production code varies.
pub(super) fn hold_instrument_constant(m: &Materialized) -> Result<()> {
    if m.baseline_dir == m.candidate_dir {
        return Ok(());
    }
    for (label, rel) in BENCH_DIRS {
        let src = m.candidate_dir.join(rel);
        let dst = m.baseline_dir.join(rel);
        std::fs::create_dir_all(&dst)?;
        copy_dir(&src, &dst).with_context(|| format!("hold {label} instrument constant"))?;
    }
    Ok(())
}
