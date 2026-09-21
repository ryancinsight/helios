//! Local replicated benchmark instrument (ADR 0003, revision 2026-09-21).
//!
//! Benchmarks are a local instrument: shared-runner wall-clock timings are
//! noise, not evidence, so CI only smoke-runs benches and this runner is the
//! sole producer of timing evidence. It executes the paired `A B B A` /
//! `B A A B` schedule on one controlled host, holding the candidate
//! benchmark sources constant across both revisions, and classifies the four
//! retained comparison roots through the Atlas-owned
//! `tools/criterion-regression` gate.
//!
//! Full default-measurement legs cannot fit a bounded suite: ten benchmark
//! IDs at Criterion defaults cost ~80 s per target per leg, which is exactly
//! why the deleted hosted job died at its 60-minute cap. This runner sizes
//! the instrument to a derived suite bound instead (benchmark time-budget
//! doctrine): abbreviated `--warm-up-time` / `--measurement-time` values
//! keep every ID and every leg, while the suite bound defaults to 1500 s --
//! the measured 1101 s of a same-revision calibration run on the reference
//! host (floored by ~15 s single iterations in `projection_throughput` that
//! no sampling abbreviation can shrink) with headroom for slower hosts.
//! Shorter measurement widens intervals and costs power against small
//! regressions; the classifier still fails closed on incomplete evidence.

use anyhow::{bail, Context, Result};
use clap::Args;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

/// The benchmark universe: `package:bench` pairs, the single source shared
/// with the CI smoke job so a binary cannot enter local measurement unsmoked.
pub const BENCHMARK_TARGETS: &[(&str, &str)] = &[
    ("helios-analysis", "dvh_queries"),
    ("helios-gpu", "projection_throughput"),
    ("helios-gpu", "transmission_throughput"),
    ("helios-solver", "scatter_superposition"),
];

/// Bench source trees held constant across both revisions.
const BENCH_DIRS: &[(&str, &str)] = &[
    ("helios-analysis", "crates/helios-analysis/benches"),
    ("helios-gpu", "crates/helios-gpu/benches"),
    ("helios-solver", "crates/helios-solver/benches"),
];

/// One replication pair: which revision runs first determines the comparison
/// direction. Block order and pair order reproduce `A B B A B A A B`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReplicationPair {
    /// `first` or `second` replication block.
    pub block: &'static str,
    /// `baseline-first` or `candidate-first` within the block.
    pub pair: &'static str,
    /// True when the baseline revision runs first in this pair.
    pub baseline_first: bool,
}

/// The committed schedule: two phase-reversed pairs per block, two blocks.
pub const SCHEDULE: &[ReplicationPair] = &[
    ReplicationPair {
        block: "first",
        pair: "baseline-first",
        baseline_first: true,
    },
    ReplicationPair {
        block: "first",
        pair: "candidate-first",
        baseline_first: false,
    },
    ReplicationPair {
        block: "second",
        pair: "candidate-first",
        baseline_first: false,
    },
    ReplicationPair {
        block: "second",
        pair: "baseline-first",
        baseline_first: true,
    },
];

/// Revision order the schedule visits: `A B B A B A A B`, balanced like the
/// deleted hosted gate (each revision occupies positions summing to 18 with
/// squared sum 102, balancing constant, linear, and quadratic period terms).
pub fn revision_sequence() -> [char; 8] {
    let mut out = ['A'; 8];
    let mut i = 0;
    for pair in SCHEDULE {
        let (first, second) = if pair.baseline_first {
            ('A', 'B')
        } else {
            ('B', 'A')
        };
        out[i] = first;
        out[i + 1] = second;
        i += 2;
    }
    out
}

#[derive(Args, Debug)]
pub struct BenchReplicatedArgs {
    /// Baseline revision (anything `git worktree add --detach` accepts).
    #[arg(long, default_value = "HEAD~1")]
    pub baseline: String,
    /// Candidate revision.
    #[arg(long, default_value = "HEAD")]
    pub candidate: String,
    /// Atlas checkout carrying `tools/criterion-regression`.
    /// Defaults to the `atlas` sibling of the helios checkout.
    #[arg(long)]
    pub atlas: Option<PathBuf>,
    /// Which replication block(s) to run. A verdict requires `both`.
    #[arg(long, value_parser = ["first", "second", "both"], default_value = "both")]
    pub pairs: String,
    /// Criterion warm-up time per benchmark ID, in seconds. Sized with
    /// `--measurement-time` so the full eight-leg schedule fits the
    /// `--suite-budget-secs` bound at the committed ten-ID universe.
    #[arg(long, default_value_t = 0.5)]
    pub warm_up_time: f64,
    /// Criterion measurement time per benchmark ID, in seconds (see
    /// `--warm-up-time` for the sizing argument).
    #[arg(long, default_value_t = 1.0)]
    pub measurement_time: f64,
    /// Per-invocation backstop for one `cargo bench` call, in seconds.
    /// Checked against measured elapsed time after the call returns; a
    /// breach fails the run loudly instead of reporting a verdict.
    #[arg(long, default_value_t = 300)]
    pub leg_timeout_secs: u64,
    /// Committed suite bound for total measurement wall-clock, in seconds.
    /// Derived, not chosen: a same-revision calibration run of the full
    /// eight-leg schedule on the reference host measured 1101 s, floored by
    /// the slowest single iterations (`projection_throughput` needs ~15 s
    /// per iteration, and Criterion collects no fewer than one iteration
    /// per sample no matter how short `--measurement-time` is). 1500 s is
    /// that measurement with headroom for slower hosts. Breaching it fails
    /// loudly: the instrument is undersized, never the production code, so
    /// the run reports the breach instead of a verdict.
    #[arg(long, default_value_t = 1500)]
    pub suite_budget_secs: u64,
    /// Reuse an existing baseline checkout instead of creating an ephemeral
    /// worktree. Must be given together with `--candidate-dir`.
    /// Same-dir/same-rev is a pipeline self-test only, never evidence.
    #[arg(long)]
    pub baseline_dir: Option<PathBuf>,
    /// Reuse an existing candidate checkout (see `--baseline-dir`).
    #[arg(long)]
    pub candidate_dir: Option<PathBuf>,
    /// Print the schedule and exit without measuring.
    #[arg(long, default_value_t = false)]
    pub plan: bool,
}

/// Run `cmd` to completion, returning its wall-clock. A non-zero exit fails
/// loudly; the caller checks the elapsed time against its budget.
fn run(cmd: &mut Command, what: &str) -> Result<Duration> {
    let start = Instant::now();
    let status = cmd.status().with_context(|| format!("spawn {what}"))?;
    let elapsed = start.elapsed();
    if !status.success() {
        bail!("{what} failed with {status} after {elapsed:?}");
    }
    Ok(elapsed)
}

fn check_budget(elapsed: Duration, budget_secs: u64, what: &str) -> Result<()> {
    if elapsed.as_secs() > budget_secs {
        bail!(
            "{what} breached its {budget_secs}s budget at {}s. The instrument \
             is undersized for this host -- no verdict is reported.",
            elapsed.as_secs()
        );
    }
    Ok(())
}

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

struct Materialized {
    baseline_dir: PathBuf,
    candidate_dir: PathBuf,
    /// Owned temp root (ephemeral worktrees); `None` for explicit dirs,
    /// whose lifetime the operator owns.
    _temp_root: Option<TempRoot>,
    run_root: PathBuf,
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

fn materialize(args: &BenchReplicatedArgs, workspace_root: &Path) -> Result<Materialized> {
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
fn hold_instrument_constant(m: &Materialized) -> Result<()> {
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

fn bench_invocation(
    manifest_dir: &Path,
    target_dir: &Path,
    package: &str,
    bench: &str,
    extra: &[String],
    timeout_secs: u64,
) -> Result<()> {
    let mut cmd = Command::new("cargo");
    // No `--locked`: the committed lock drifts with moving provider heads
    // and the repo's own tests/docs jobs build unlocked; lock integrity is
    // the separate lockfile-guard job's contract, not this instrument's.
    // `CARGO_TARGET_DIR` in the environment is load-bearing beyond the
    // `--target-dir` flag: Criterion child processes resolve their report
    // root from the environment, and without it they write into the shared
    // stack cache beside other repos' benchmark data. Both point at the
    // run-private target dir.
    cmd.env("CARGO_TARGET_DIR", target_dir);
    cmd.arg("bench")
        .arg("--manifest-path")
        .arg(manifest_dir.join("Cargo.toml"))
        .arg("--target-dir")
        .arg(target_dir)
        .arg("--all-features")
        .arg("--package")
        .arg(package)
        .arg("--bench")
        .arg(bench)
        .arg("--")
        .args(extra);
    let elapsed = run(
        &mut cmd,
        &format!("cargo bench {package}/{bench} {extra:?}"),
    )?;
    check_budget(
        elapsed,
        timeout_secs,
        &format!("cargo bench {package}/{bench}"),
    )
}

fn smoke_candidate(m: &Materialized, target_dir: &Path) -> Result<()> {
    // The smoke gate precedes every leg: new binaries cannot enter
    // regression measurement unsmoked. The 120 s per-target bound derives
    // from a measured 71 s single-iteration sweep of the slowest target
    // (`projection_throughput`, 4 IDs) on the calibration host, with
    // headroom for slower hosts.
    for (package, bench) in BENCHMARK_TARGETS {
        bench_invocation(
            &m.candidate_dir,
            target_dir,
            package,
            bench,
            &["--test".to_string()],
            120,
        )?;
    }
    Ok(())
}

fn full_leg(
    manifest_dir: &Path,
    target_dir: &Path,
    args: &BenchReplicatedArgs,
    baseline_name: &str,
    confidence: Option<&str>,
) -> Result<()> {
    for (package, bench) in BENCHMARK_TARGETS {
        let mut extra = vec![
            "--warm-up-time".to_string(),
            args.warm_up_time.to_string(),
            "--measurement-time".to_string(),
            args.measurement_time.to_string(),
        ];
        if let Some(level) = confidence {
            extra.push("--baseline".to_string());
            extra.push(baseline_name.to_string());
            extra.push("--confidence-level".to_string());
            extra.push(level.to_string());
        } else {
            extra.push("--save-baseline".to_string());
            extra.push(baseline_name.to_string());
        }
        bench_invocation(
            manifest_dir,
            target_dir,
            package,
            bench,
            &extra,
            args.leg_timeout_secs,
        )?;
    }
    Ok(())
}

/// Derive the family-wise confidence from the first leg's criterion root.
/// No `--locked`: under the stack checkout the shared `[patch]` overlay
/// forces re-resolution of any manifest it touches, while the Atlas tool
/// lockfiles stay overlay-clean by policy; the gate builds unlocked.
fn required_confidence(
    atlas_tool: &Path,
    target_dir: &Path,
    baseline_name: &str,
) -> Result<String> {
    let output = Command::new("cargo")
        .arg("run")
        .arg("--manifest-path")
        .arg(atlas_tool.join("tools/criterion-regression/Cargo.toml"))
        .arg("--target-dir")
        .arg(target_dir.join("gate"))
        .arg("--")
        .arg("required-confidence")
        .arg("--criterion-root")
        .arg(target_dir.join("criterion"))
        .arg("--baseline")
        .arg(baseline_name)
        .output()
        .context("run required-confidence")?;
    if !output.status.success() {
        bail!("required-confidence failed with {}", output.status);
    }
    Ok(String::from_utf8(output.stdout)?.trim().to_string())
}

fn classify(
    atlas_tool: &Path,
    target_dir: &Path,
    reports: &Path,
    baseline_name: &str,
) -> Result<()> {
    // The gate exits 0 on a clean verdict, 1 on a verdict with failures,
    // and 2 on a tool error. Exit 1 is a verdict, not a crash: surface the
    // `replicated result:` line and keep a non-zero exit so automation
    // notices, without mislabeling evidence as a tool failure.
    let output = Command::new("cargo")
        .arg("run")
        .arg("--manifest-path")
        .arg(atlas_tool.join("tools/criterion-regression/Cargo.toml"))
        .arg("--target-dir")
        .arg(target_dir.join("gate"))
        .arg("--")
        .arg("check-replicated-counterbalanced")
        .arg("--first-baseline-first-root")
        .arg(reports.join("first").join("baseline-first"))
        .arg("--first-candidate-first-root")
        .arg(reports.join("first").join("candidate-first"))
        .arg("--second-baseline-first-root")
        .arg(reports.join("second").join("baseline-first"))
        .arg("--second-candidate-first-root")
        .arg(reports.join("second").join("candidate-first"))
        .arg("--baseline")
        .arg(baseline_name)
        .output()
        .context("run check-replicated-counterbalanced")?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    interpret_verdict(output.status.code(), &stdout, &stderr, reports)
}

/// Pure verdict interpretation, unit-tested below: exit 0 is a clean
/// verdict, exit 1 a verdict with failures, anything else a tool failure.
fn interpret_verdict(code: Option<i32>, stdout: &str, stderr: &str, reports: &Path) -> Result<()> {
    for line in stdout.lines() {
        if line.starts_with("replicated regression:") || line.starts_with("replicated result:") {
            println!("{line}");
        }
    }
    println!("reports retained under {}", reports.display());
    match code {
        Some(0) => Ok(()),
        Some(1) => bail!(
            "verdict: replicated regression(s) or universe mismatch(es); evidence retained under {}. \
             A verdict on identical revisions names host noise, never a code defect -- re-run on a \
             controlled host before attributing it to production code.",
            reports.display()
        ),
        _ => bail!(
            "check-replicated-counterbalanced failed with exit code {code:?}. stderr: {}",
            stderr.trim()
        ),
    }
}

/// Human-readable schedule plan; also the `--plan` output. The revision
/// order line derives from [`revision_sequence`], so the plan and the
/// executed schedule cannot drift apart.
pub fn plan_text(baseline: &str, candidate: &str, pairs: &str) -> String {
    use std::fmt::Write as _;
    let mut out = format!(
        "paired schedule {} (baseline {baseline} vs candidate {candidate}):\n",
        revision_sequence().iter().collect::<String>()
    );
    for pair in SCHEDULE {
        if pairs != "both" && pair.block != pairs {
            continue;
        }
        let (first, second) = if pair.baseline_first {
            ("baseline", "candidate")
        } else {
            ("candidate", "baseline")
        };
        let _ = writeln!(out, "  {} {}: {first} then {second}", pair.block, pair.pair);
    }
    out
}

pub fn run_bench_replicated(args: &BenchReplicatedArgs, workspace_root: &Path) -> Result<()> {
    print!(
        "{}",
        plan_text(&args.baseline, &args.candidate, &args.pairs)
    );
    if args.plan {
        return Ok(());
    }
    let atlas_tool = match &args.atlas {
        Some(path) => path.clone(),
        None => workspace_root.parent().map_or_else(
            || workspace_root.join("..").join("atlas"),
            |p| p.join("atlas"),
        ),
    };
    if !atlas_tool
        .join("tools/criterion-regression/Cargo.toml")
        .exists()
    {
        bail!(
            "Atlas criterion-regression gate not found under {}; pass --atlas <checkout>",
            atlas_tool.display()
        );
    }

    let m = materialize(args, workspace_root)?;
    let target_dir = m.run_root.join("target");
    let reports = m.run_root.join("reports");
    let baseline_name = "atlas-base";

    hold_instrument_constant(&m)?;
    smoke_candidate(&m, &target_dir)?;

    let selected: Vec<ReplicationPair> = SCHEDULE
        .iter()
        .copied()
        .filter(|pair| args.pairs == "both" || pair.block == args.pairs)
        .collect();

    let suite_start = Instant::now();
    let mut confidence: Option<String> = None;
    for pair in &selected {
        let pair_start = Instant::now();
        let (first_dir, second_dir) = if pair.baseline_first {
            (&m.baseline_dir, &m.candidate_dir)
        } else {
            (&m.candidate_dir, &m.baseline_dir)
        };
        let criterion_root = target_dir.join("criterion");
        if criterion_root.exists() {
            std::fs::remove_dir_all(&criterion_root)?;
        }
        full_leg(first_dir, &target_dir, args, baseline_name, None)?;
        if confidence.is_none() {
            confidence = Some(required_confidence(
                &atlas_tool,
                &target_dir,
                baseline_name,
            )?);
        }
        full_leg(
            second_dir,
            &target_dir,
            args,
            baseline_name,
            confidence.as_deref(),
        )?;
        let dest = reports.join(pair.block).join(pair.pair);
        if dest.exists() {
            std::fs::remove_dir_all(&dest)?;
        }
        let Some(parent) = dest.parent() else {
            bail!("report destination {} has no parent", dest.display());
        };
        std::fs::create_dir_all(parent)?;
        std::fs::rename(&criterion_root, &dest)
            .with_context(|| format!("retain {} {} comparison root", pair.block, pair.pair))?;
        println!(
            "pair {} {}: {}s",
            pair.block,
            pair.pair,
            pair_start.elapsed().as_secs()
        );
    }
    let elapsed_secs = suite_start.elapsed().as_secs();
    println!(
        "measurement wall-clock: {elapsed_secs}s (suite bound: {}s)",
        args.suite_budget_secs
    );
    if elapsed_secs > args.suite_budget_secs {
        bail!(
            "suite bound breached: {elapsed_secs}s > {}s. The instrument is \
             undersized for this host -- reduce --measurement-time / \
             --warm-up-time or raise --suite-budget-secs with a recorded \
             justification. No verdict is reported.",
            args.suite_budget_secs
        );
    }

    if args.pairs == "both" {
        classify(&atlas_tool, &target_dir, &reports, baseline_name)?;
        println!("reports retained under {}", reports.display());
    } else {
        println!(
            "partial schedule (--pairs {}): roots retained under {}; \
             a verdict requires --pairs both",
            args.pairs,
            reports.display()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schedule_visits_a_b_b_a_b_a_a_b() {
        assert_eq!(
            revision_sequence(),
            ['A', 'B', 'B', 'A', 'B', 'A', 'A', 'B']
        );
    }

    #[test]
    fn schedule_balances_each_revision_like_the_hosted_gate() {
        // Positions are 1-based; each revision sums to 18 with squared sum
        // 102, balancing exposure to constant, linear, and quadratic period
        // terms -- the same balance the deleted hosted gate documented.
        for want in ['A', 'B'] {
            let positions: Vec<u32> = revision_sequence()
                .iter()
                .enumerate()
                .filter(|(_, r)| **r == want)
                .map(|(i, _)| (i + 1) as u32)
                .collect();
            assert_eq!(positions.iter().sum::<u32>(), 18);
            assert_eq!(positions.iter().map(|p| p * p).sum::<u32>(), 102);
        }
    }

    #[test]
    fn plan_lists_both_blocks_by_default() {
        let text = plan_text("HEAD~1", "HEAD", "both");
        assert!(text.contains("first baseline-first"));
        assert!(text.contains("second candidate-first"));
    }

    #[test]
    fn plan_filters_to_one_block() {
        let text = plan_text("HEAD~1", "HEAD", "first");
        assert!(text.contains("first baseline-first"));
        assert!(!text.contains("second"));
    }

    #[test]
    fn every_target_names_a_nonempty_package_and_bench() {
        // The CI smoke job and this runner share BENCHMARK_TARGETS through
        // this file; the smoke list is asserted identical by construction.
        assert_eq!(BENCHMARK_TARGETS.len(), 4);
        for (package, bench) in BENCHMARK_TARGETS {
            assert!(!package.is_empty() && !bench.is_empty());
        }
    }

    #[test]
    fn verdict_clean_is_ok() {
        let reports = Path::new("reports");
        let out = "replicated result: 0 regression(s), 0 replication-universe mismatch(es)\n";
        assert!(interpret_verdict(Some(0), out, "", reports).is_ok());
    }

    #[test]
    fn verdict_with_failures_is_a_verdict_not_a_crash() {
        let reports = Path::new("reports");
        let out = "replicated regression: foo first +1.00%/+2.00%; second +3.00%/+4.00%\n\
                   replicated result: 1 regression(s), 0 replication-universe mismatch(es)\n";
        let err = interpret_verdict(Some(1), out, "", reports)
            .expect_err("exit 1 carries a verdict, never Ok");
        assert!(err.to_string().starts_with("verdict:"), "{err}");
    }

    #[test]
    fn verdict_tool_error_is_a_failure() {
        let reports = Path::new("reports");
        let err = interpret_verdict(Some(2), "", "boom", reports)
            .expect_err("exit 2 is a tool error, never Ok");
        assert!(err.to_string().contains("failed"), "{err}");
    }
}
