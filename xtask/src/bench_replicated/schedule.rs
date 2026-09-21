//! Counterbalanced revision schedule and its rendered plan.

/// One replication pair: which revision runs first determines the comparison
/// direction. Block order and pair order reproduce `A B B A B A A B`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct ReplicationPair {
    /// `first` or `second` replication block.
    pub block: &'static str,
    /// `baseline-first` or `candidate-first` within the block.
    pub pair: &'static str,
    /// True when the baseline revision runs first in this pair.
    pub baseline_first: bool,
}

/// The committed schedule: two phase-reversed pairs per block, two blocks.
pub(super) const SCHEDULE: &[ReplicationPair] = &[
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
pub(super) fn revision_sequence() -> [char; 8] {
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

/// Human-readable schedule plan; also the `--plan` output. The revision
/// order line derives from [`revision_sequence`], so the plan and the
/// executed schedule cannot drift apart.
pub(super) fn plan_text(baseline: &str, candidate: &str, pairs: &str) -> String {
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
}
