//! Bounded benchmark execution and replicated-verdict classification.

use super::{checkout::Materialized, BenchReplicatedArgs};
use anyhow::{bail, Context, Result};
use std::path::Path;
use std::process::Command;
use std::time::{Duration, Instant};

/// The benchmark universe: `package:bench` pairs, the single source shared
/// with the CI smoke job so a binary cannot enter local measurement unsmoked.
const BENCHMARK_TARGETS: &[(&str, &str)] = &[
    ("helios-analysis", "dvh_queries"),
    ("helios-gpu", "projection_throughput"),
    ("helios-gpu", "transmission_throughput"),
    ("helios-imaging", "ramp_filter"),
    ("helios-solver", "scatter_superposition"),
];

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

pub(super) fn smoke_candidate(m: &Materialized, target_dir: &Path) -> Result<()> {
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

pub(super) fn full_leg(
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
pub(super) fn required_confidence(
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

pub(super) fn classify(
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_target_names_a_nonempty_package_and_bench() {
        // The CI smoke job and this runner share BENCHMARK_TARGETS through
        // this file; the smoke list is asserted identical by construction.
        assert_eq!(BENCHMARK_TARGETS.len(), 5);
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
