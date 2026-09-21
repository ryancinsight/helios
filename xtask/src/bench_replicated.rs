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
use std::time::Instant;

mod checkout;
mod measurement;
mod schedule;

use self::checkout::{hold_instrument_constant, materialize};
use self::measurement::{classify, full_leg, required_confidence, smoke_candidate};
use self::schedule::{plan_text, ReplicationPair, SCHEDULE};

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
