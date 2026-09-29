    use super::*;
    use eunomia::assert_relative_eq;

    fn dose(value: f64) -> AbsorbedDose<f64> {
        AbsorbedDose::from_base(value)
    }

    fn dimensionless(value: f64) -> Dimensionless<f64> {
        Dimensionless::from_base(value)
    }

    // Independent scalar gEUD oracle for the coeus-tape gEUD — deliberately a
    // separate implementation from both the tape and helios-analysis's
    // `generalized_eud` (a differential test must not check code against itself).
    fn geud_ref(doses: &[f64], a: f64) -> f64 {
        let n = doses.len() as f64;
        (doses.iter().map(|&d| d.powf(a)).sum::<f64>() / n).powf(1.0 / a)
    }

    /// 3×2 influence with distinct entries; the differential oracle is the exact
    /// hand gradient Aᵀ(A·x − d) from the projected-gradient solver.
    fn influence() -> DoseInfluence<f64> {
        DoseInfluence::from_rows(3, 2, vec![1.0, 2.0, 0.5, -1.0, 3.0, 4.0])
            .expect("synthetic dose-influence matrix is well-formed")
    }

    // All-positive 3×2 influence so `A·x > 0` for positive `x` (gEUD with a
    // non-integer `a` needs positive dose).
    fn positive_influence() -> DoseInfluence<f64> {
        DoseInfluence::from_rows(3, 2, vec![1.0, 0.5, 2.0, 1.0, 0.5, 1.5])
            .expect("synthetic dose-influence matrix is well-formed")
    }

    #[test]
    fn tape_geud_value_matches_the_analytic_geud() {
        // UpperLimit with reference 0 makes the objective L = gEUD², so √L is the
        // tape's gEUD — cross-check it against the independent geud_ref oracle.
        let inf = positive_influence();
        let x = [0.8, 0.6];
        let a = 2.5;
        let penalty = EudPenalty {
            a: dimensionless(a),
            reference: dose(0.0),
            kind: EudKind::UpperLimit,
            weight: 1.0,
        };
        let (value, _) = eud_objective_gradient_autodiff(&inf, &x, &penalty)
            .expect("planning fixture round-trips");
        let analytic = geud_ref(&inf.apply(&x), a);
        assert_relative_eq!(value.sqrt(), analytic, max_relative = 1e-10);
    }

    #[test]
    fn eud_gradient_matches_central_finite_difference() {
        // Active upper hinge: pin the tape gradient against a central finite
        // difference of the objective (the differential oracle over the whole
        // gEUD-plus-penalty tape).
        let inf = positive_influence();
        let x = [0.8, 0.6];
        let a = 3.0;
        let geud = geud_ref(&inf.apply(&x), a);
        let penalty = EudPenalty {
            a: dimensionless(a),
            reference: dose(0.5 * geud), // strictly below gEUD ⇒ hinge active
            kind: EudKind::UpperLimit,
            weight: 2.0,
        };
        let value_at = |x: &[f64]| {
            eud_objective_gradient_autodiff(&inf, x, &penalty)
                .expect("finite sample points keep the objective solvable")
                .0
        };
        let (_, grad) = eud_objective_gradient_autodiff(&inf, &x, &penalty)
            .expect("planning fixture round-trips");
        let h = 1e-6;
        for k in 0..x.len() {
            let (mut xp, mut xm) = (x, x);
            xp[k] += h;
            xm[k] -= h;
            let fd = (value_at(&xp) - value_at(&xm)) / (2.0 * h);
            assert_relative_eq!(grad[k], fd, max_relative = 1e-5, epsilon = 1e-7);
        }
    }

    #[test]
    fn eud_gradient_is_zero_within_the_limit() {
        // gEUD strictly below the upper reference ⇒ hinge inactive ⇒ L = 0, ∇ = 0.
        let inf = positive_influence();
        let x = [0.8, 0.6];
        let a = 3.0;
        let geud = geud_ref(&inf.apply(&x), a);
        let penalty = EudPenalty {
            a: dimensionless(a),
            reference: dose(geud * 2.0), // well above ⇒ no violation
            kind: EudKind::UpperLimit,
            weight: 5.0,
        };
        let (value, grad) = eud_objective_gradient_autodiff(&inf, &x, &penalty)
            .expect("planning fixture round-trips");
        assert_relative_eq!(value, 0.0, epsilon = 1e-12);
        for g in grad {
            assert_relative_eq!(g, 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn eud_zero_a_and_shape_mismatch_are_typed_errors() {
        let inf = positive_influence();
        let bad_a = EudPenalty {
            a: dimensionless(0.0),
            reference: dose(1.0),
            kind: EudKind::UpperLimit,
            weight: 1.0,
        };
        assert_eq!(
            eud_objective_gradient_autodiff(&inf, &[0.8, 0.6], &bad_a),
            Err(HeliosError::InvalidDomainValue {
                field: "eud_objective_gradient_autodiff::a",
                value: 0.0,
                reason: "gEUD volume parameter must be finite and non-zero",
            })
        );
        let ok = EudPenalty {
            a: dimensionless(2.0),
            reference: dose(1.0),
            kind: EudKind::UpperLimit,
            weight: 1.0,
        };
        assert_eq!(
            eud_objective_gradient_autodiff(&inf, &[0.8], &ok),
            Err(HeliosError::InvalidDomainValue {
                field: "eud_objective_gradient_autodiff::x",
                value: 1.0,
                reason: "weight count must equal the beamlet count",
            })
        );
    }

    #[test]
    fn autodiff_gradient_matches_the_exact_hand_gradient() {
        let inf = influence();
        let x = [0.7, -0.3];
        let d = [1.0, 0.5, 2.0];

        // Exact: Aᵀ(A·x − d).
        let ax = inf.apply(&x);
        let residual: Vec<f64> = ax.iter().zip(&d).map(|(&a, &b)| a - b).collect();
        let exact = inf.transpose_apply(&residual);

        let auto = objective_gradient_autodiff(&inf, &x, &d).expect("gradient");
        assert_eq!(auto.len(), exact.len());
        for (i, (&g, &e)) in auto.iter().zip(&exact).enumerate() {
            // Same values through a different summation route: bound at 1e-12
            // relative (f64; tiny fixed dims, reduction depth ≤ 3).
            assert_relative_eq!(g, e, max_relative = 1e-12, epsilon = 1e-14);
            let _ = i;
        }
    }

    #[test]
    fn gradient_is_zero_at_the_least_squares_optimum() {
        // Identity influence: optimum x* = d exactly → ∇ = 0.
        let inf = DoseInfluence::from_rows(2, 2, vec![1.0, 0.0, 0.0, 1.0])
            .expect("planning fixture round-trips");
        let d = [2.0, -1.5];
        let grad = objective_gradient_autodiff(&inf, &d, &d).expect("gradient");
        for &g in &grad {
            assert_relative_eq!(g, 0.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn dvh_gradient_matches_the_hand_subgradient() {
        // ∇L = −2·w_u·Aᵀ relu(f − Ax) + 2·w_o·Aᵀ relu(Ax − c). Choose x so both
        // hinges are strictly active somewhere (no kink ambiguity).
        let inf = influence();
        let x = [0.4, -0.2];
        let floor = [dose(1.0), dose(0.2), dose(0.5)];
        let ceiling = [dose(1.5), dose(0.4), dose(0.9)];
        let penalty = DvhPenalty {
            floor: &floor,
            ceiling: &ceiling,
            weight_under: 2.0,
            weight_over: 3.0,
        };
        let ax = inf.apply(&x);
        let floor_values: Vec<f64> = floor.iter().map(|value| *value.as_base()).collect();
        let ceiling_values: Vec<f64> = ceiling.iter().map(|value| *value.as_base()).collect();
        let under: Vec<f64> = ax
            .iter()
            .zip(&floor_values)
            .map(|(&d, &f)| (f - d).max(0.0))
            .collect();
        let over: Vec<f64> = ax
            .iter()
            .zip(&ceiling_values)
            .map(|(&d, &c)| (d - c).max(0.0))
            .collect();
        let gu = inf.transpose_apply(&under);
        let go = inf.transpose_apply(&over);
        let hand: Vec<f64> = gu
            .iter()
            .zip(&go)
            .map(|(&u, &o)| -2.0 * 2.0 * u + 2.0 * 3.0 * o)
            .collect();

        let (value, auto) = dvh_objective_gradient_autodiff(&inf, &x, &penalty).expect("grad");
        // Objective value cross-check.
        let expected_value = 2.0 * under.iter().map(|u| u * u).sum::<f64>()
            + 3.0 * over.iter().map(|o| o * o).sum::<f64>();
        assert_relative_eq!(value, expected_value, max_relative = 1e-12);
        for (&g, &h) in auto.iter().zip(&hand) {
            assert_relative_eq!(g, h, max_relative = 1e-12, epsilon = 1e-14);
        }
    }

    #[test]
    fn gradient_is_zero_inside_the_penalty_band() {
        // Dose strictly between floor and ceiling everywhere → both hinges
        // inactive → L = 0 and ∇L = 0 (the band is free).
        let inf = DoseInfluence::from_rows(2, 2, vec![1.0, 0.0, 0.0, 1.0])
            .expect("planning fixture round-trips");
        let x = [1.0, 1.0]; // dose = [1, 1]
        let penalty = DvhPenalty {
            floor: &[dose(0.5), dose(0.5)],
            ceiling: &[dose(2.0), dose(2.0)],
            weight_under: 1.0,
            weight_over: 1.0,
        };
        let (value, grad) = dvh_objective_gradient_autodiff(&inf, &x, &penalty).expect("grad");
        assert_relative_eq!(value, 0.0, epsilon = 1e-14);
        for &g in &grad {
            assert_relative_eq!(g, 0.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn dvh_optimizer_meets_target_floor_and_oar_ceiling() {
        // Voxel 0 is the target (floor 1.0), voxel 1 the OAR (ceiling 0.3).
        // Beamlet 0 doses both (target 1.0, OAR 0.5/unit); beamlet 1 doses only
        // the target. The optimum uses beamlet 1 (OAR-sparing) to reach the floor
        // while keeping the OAR under its ceiling.
        let inf = DoseInfluence::from_rows(2, 2, vec![1.0, 1.0, 0.5, 0.0])
            .expect("planning fixture round-trips");
        let penalty = DvhPenalty {
            floor: &[dose(1.0), dose(0.0)],
            ceiling: &[dose(10.0), dose(0.3)],
            weight_under: 1.0,
            weight_over: 10.0,
        };
        let x = optimize_beam_weights_dvh(&inf, &penalty, 500, 0.05).expect("optimize");
        let dose = inf.apply(&x);
        assert!(x.iter().all(|&w| w >= 0.0), "weights must be non-negative");
        assert!(dose[0] > 0.95, "target dose {} below floor", dose[0]);
        assert!(dose[1] < 0.35, "OAR dose {} above ceiling", dose[1]);
    }

    #[test]
    fn shape_mismatches_are_typed_errors() {
        let inf = influence();
        assert!(objective_gradient_autodiff(&inf, &[1.0], &[1.0, 1.0, 1.0]).is_err());
        assert!(objective_gradient_autodiff(&inf, &[1.0, 2.0], &[1.0]).is_err());
    }
