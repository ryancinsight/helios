//! End-to-end helical TomoTherapy + MVCT workflow validation (H-041).
//!
//! One shared attenuation volume `μ` drives BOTH branches of the platform, proving
//! the layers compose across their seams:
//!
//! - **Imaging / IGRT:** `μ` → parallel-beam Radon → FBP reconstruction (recovers
//!   `μ` in a water ROI) → rigid registration recovers a known couch shift.
//! - **Therapy:** helical MLC delivery → per-leaf divergent-fan terma deposition
//!   into `μ` → collapsed-cone scatter → dose → DVH + 3%/2 mm gamma self-consistency.
//!
//! Every assertion is an analytical / self-consistency oracle; no external engine
//! or licensed dataset is involved (those gates are environment-blocked, see
//! gap_audit G-16/G-18).

use aequitas::systems::si::{
    quantities::{AbsorbedDose, AreaPerMass, Length},
    units::{Millimeter, SquareCentimeterPerGram},
};
use helios_analysis::{gamma_index_3d, gamma_pass_rate, roi_statistics, spherical_mask, Dvh};
use helios_domain::{HelicalDelivery, LeafOpenTimeSinogram, MlcModel, Volume, VoxelGrid};
use helios_imaging::{filtered_back_projection, parallel_beam_radon, register_translation};
use helios_math::Point3;
use helios_simulation::{
    accumulate_delivered_dose, accumulate_delivered_dose_anisotropic, simulate_helical_delivery,
    BeamGeometry, CollapsedCone, SpectralComponent,
};
use helios_solver::{attenuation_map, scatter_superposition, symmetric_deposition_kernel};
use hyperion::coefficient::MassAttenuation;

const NX: usize = 31;
const NZ: usize = 5;
const SPACING: f64 = 2.0;
const CENTRE_MM: f64 = (NX as f64 - 1.0) * SPACING / 2.0; // 30 mm

fn water_mass_attenuation() -> MassAttenuation<f64> {
    MassAttenuation::new(AreaPerMass::from_unit::<SquareCentimeterPerGram>(0.06))
        .expect("valid water mass attenuation")
}

// CT phantom (HU): air outside a 24 mm water cylinder, with an 8 mm bone insert
// (HU 800) at the centre. Uniform along z (cylinder axis).
fn ct_phantom() -> Volume<f64> {
    let grid =
        VoxelGrid::axis_aligned([NX, NX, NZ], [SPACING; 3], Point3::new(0.0, 0.0, 0.0)).unwrap();
    Volume::from_shape_fn(grid, |idx| {
        let dx = idx[0] as f64 * SPACING - CENTRE_MM;
        let dy = idx[1] as f64 * SPACING - CENTRE_MM;
        let r = (dx * dx + dy * dy).sqrt();
        if r <= 8.0 {
            800.0 // bone insert
        } else if r <= 24.0 {
            0.0 // water
        } else {
            -1000.0 // air
        }
    })
}

#[test]
fn shared_mu_drives_imaging_and_delivery_end_to_end() {
    // ── Shared physics: CT → μ (μ/ρ = 0.06 cm²/g, water 1.0 g/cm³). ──
    let ct = ct_phantom();
    let mu =
        attenuation_map(&ct, water_mass_attenuation(), 1.0).expect("fixture calibration is finite");
    // Sanity: water μ = 0.06·1.0; bone μ = 0.06·1.8 (HU 800 → density 1.8); air ≈ 0.
    assert!((mu.get(15, 15, 2).unwrap() - 0.108).abs() < 1e-9); // centre = bone
    assert!(mu.get(0, 0, 2).unwrap().abs() < 1e-6); // corner = air

    // ── Imaging branch: Radon → FBP; recover μ in a water ROI. ──
    let angles: Vec<f64> = (0..90)
        .map(|a| a as f64 * std::f64::consts::PI / 90.0)
        .collect();
    let n_off = 61;
    let offsets: Vec<f64> = (0..n_off)
        .map(|j| -30.0 + j as f64 * 60.0 / (n_off - 1) as f64)
        .collect();
    let sino = parallel_beam_radon(&mu, &angles, &offsets, 400.0, 0.5);
    let recon_grid =
        VoxelGrid::axis_aligned([NX, NX, 1], [SPACING; 3], Point3::new(0.0, 0.0, 0.0)).unwrap();
    let recon = filtered_back_projection(&sino, &recon_grid);
    // Water ROI (well inside the cylinder, outside the insert; r ≈ 10–16 mm).
    let water = roi_statistics(&recon, [20, 14, 0], [23, 17, 1]);
    assert!(
        (water.mean - 0.06).abs() / 0.06 < 0.2,
        "recon water μ {} not within 20% of 0.06",
        water.mean
    );

    // ── IGRT: register a known couch shift on the (textured) μ centre slice. ──
    let slice = |v: &Volume<f64>, k: usize| {
        Volume::from_shape_fn(recon_grid, move |idx| v.get(idx[0], idx[1], k).unwrap())
    };
    let fixed = slice(&mu, 2);
    // moving = μ centre slice shifted by (2, −1): moving(v) = fixed(v − s).
    let moving = Volume::from_shape_fn(recon_grid, |idx| {
        let (fi, fj) = (idx[0] as isize - 2, idx[1] as isize + 1);
        if fi >= 0 && fj >= 0 && (fi as usize) < NX && (fj as usize) < NX {
            fixed.get(fi as usize, fj as usize, 0).unwrap()
        } else {
            0.0
        }
    });
    assert_eq!(register_translation(&fixed, &moving, [3, 3, 0]), [2, -1, 0]);

    // ── Therapy branch: helical MLC delivery → dose. ──
    let leaves = 9;
    let projections = 20;
    let delivery = HelicalDelivery::new(51, 20.0, 0.3, 10.0, 0.0, 0.0).unwrap();
    let lot =
        LeafOpenTimeSinogram::from_fractions(projections, leaves, vec![1.0; projections * leaves])
            .unwrap();
    let mlc = MlcModel::new(0.01, 0.05).unwrap();
    let frames = simulate_helical_delivery(&delivery, &lot, &mlc);
    assert_eq!(frames.len(), projections);

    let terma = accumulate_delivered_dose(
        &frames,
        &mu,
        BeamGeometry::PointSource {
            source_axis_mm: 850.0,
        },
        4.0,
        0.5,
    )
    .expect("attenuation map satisfies Hyperion's transport contract");
    let kernel = symmetric_deposition_kernel(0.5_f64, 0.2, 1);
    let dose = scatter_superposition(&terma, &kernel, &kernel, &kernel);

    // The delivery deposited energy into the patient, and dose is non-negative.
    assert!(dose.sum() > 0.0, "delivery deposited no dose");
    for k in 0..NZ {
        for j in 0..NX {
            for i in 0..NX {
                assert!(dose.get(i, j, k).unwrap() >= 0.0);
            }
        }
    }

    // ── Analysis: DVH + 3%/2 mm gamma self-consistency. ──
    let dvh = Dvh::from_volume(&dose);
    assert!(
        *dvh.mean().as_base() > 0.0,
        "DVH mean dose must be positive"
    );
    // Gamma of the dose against itself is identically 0 → 100% pass at 3%/2 mm.
    let peak = dose.mean_top_dose();
    let gamma = gamma_index_3d(
        &dose,
        &dose,
        0.03,
        Length::from_unit::<Millimeter>(2.0),
        AbsorbedDose::from_base(peak),
        Length::from_unit::<Millimeter>(6.0),
    )
    .unwrap();
    let pass = gamma_pass_rate(&gamma, &dose, AbsorbedDose::from_base(0.1 * peak));
    assert!(
        (pass - 1.0).abs() < 1e-9,
        "self-gamma pass rate {pass} must be 100%"
    );
}

#[test]
fn beam_following_poly_energetic_dose_end_to_end() {
    // The therapy branch through the H-020i/H-020j dose model: helical delivery →
    // per-frame terma → beam-following, poly-energetic (beam-hardened) collapsed
    // cone → dose. Oracles are analytical / self-consistent (no external engine).
    let ct = ct_phantom();
    let mu =
        attenuation_map(&ct, water_mass_attenuation(), 1.0).expect("fixture calibration is finite");

    let leaves = 9;
    let projections = 20;
    let delivery = HelicalDelivery::new(51, 20.0, 0.3, 10.0, 0.0, 0.0).unwrap();
    let lot =
        LeafOpenTimeSinogram::from_fractions(projections, leaves, vec![1.0; projections * leaves])
            .unwrap();
    let mlc = MlcModel::new(0.01, 0.05).unwrap();
    let frames = simulate_helical_delivery(&delivery, &lot, &mlc);
    let geom = BeamGeometry::PointSource {
        source_axis_mm: 850.0,
    };

    // Two-component (soft/hard) forward-peaked spectrum, re-oriented per gantry angle.
    let voxel_cm = SPACING / 10.0; // 0.2 cm
    let spectrum = [
        SpectralComponent {
            range_up_cm: 0.10,
            range_down_cm: 0.35,
            weight: 0.7,
        },
        SpectralComponent {
            range_up_cm: 0.05,
            range_down_cm: 0.90,
            weight: 0.3,
        },
    ];
    let cone = CollapsedCone::poly_forward_peaked(&spectrum, 0.5, voxel_cm, 2, 3, 1);
    let dose = accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 4.0, 0.5, &cone)
        .expect("attenuation map satisfies Hyperion's transport contract");

    // Non-negativity everywhere; the delivery deposited energy.
    assert!(dose.sum() > 0.0, "delivery deposited no dose");
    for k in 0..NZ {
        for j in 0..NX {
            for i in 0..NX {
                assert!(dose.get(i, j, k).unwrap() >= 0.0, "negative dose");
            }
        }
    }

    // Energy conservation: the Σ=1 collapsed cone can only redistribute or truncate
    // energy at the boundary, never create it — so the scattered dose total is ≤ the
    // deposited terma total, and (delivery concentrated near centre) retains most of
    // it. `accumulate_delivered_dose` with the same frames/geometry yields the terma.
    let terma = accumulate_delivered_dose(&frames, &mu, geom, 4.0, 0.5)
        .expect("attenuation map satisfies Hyperion's transport contract");
    let (dsum, tsum) = (dose.sum(), terma.sum());
    assert!(
        dsum <= tsum * (1.0 + 1e-9),
        "scattered dose {dsum} exceeds deposited terma {tsum} (energy created)"
    );
    assert!(
        dsum > tsum * 0.85,
        "scattered dose {dsum} lost too much of terma {tsum} to boundary truncation"
    );

    // DVH + 3%/2 mm self-gamma self-consistency (dose vs itself → 100% pass).
    let dvh = Dvh::from_volume(&dose);
    assert!(
        *dvh.mean().as_base() > 0.0,
        "DVH mean dose must be positive"
    );
    let peak = dose.mean_top_dose();
    let gamma = gamma_index_3d(
        &dose,
        &dose,
        0.03,
        Length::from_unit::<Millimeter>(2.0),
        AbsorbedDose::from_base(peak),
        Length::from_unit::<Millimeter>(6.0),
    )
    .unwrap();
    let pass = gamma_pass_rate(&gamma, &dose, AbsorbedDose::from_base(0.1 * peak));
    assert!(
        (pass - 1.0).abs() < 1e-9,
        "self-gamma pass rate {pass} must be 100%"
    );
}

#[test]
fn per_structure_plan_evaluation_over_delivered_dose() {
    // Run the per-structure plan-evaluation surface over real delivered dose:
    // masked DVH → gEUD → TCP (target) / NTCP (OAR). Oracles are clinical-
    // plausibility + probability well-formedness (no external engine).
    let ct = ct_phantom();
    let mu =
        attenuation_map(&ct, water_mass_attenuation(), 1.0).expect("fixture calibration is finite");
    let grid = *mu.grid();

    let (leaves, projections) = (9, 20);
    let delivery = HelicalDelivery::new(51, 20.0, 0.3, 10.0, 0.0, 0.0).unwrap();
    let lot =
        LeafOpenTimeSinogram::from_fractions(projections, leaves, vec![1.0; projections * leaves])
            .unwrap();
    let mlc = MlcModel::new(0.01, 0.05).unwrap();
    let frames = simulate_helical_delivery(&delivery, &lot, &mlc);
    let geom = BeamGeometry::PointSource {
        source_axis_mm: 850.0,
    };
    let voxel_cm = SPACING / 10.0;
    let cone = CollapsedCone::forward_peaked(0.1, 0.6, 0.5, voxel_cm, 2, 3, 1);
    let dose = accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 4.0, 0.5, &cone)
        .expect("attenuation map satisfies Hyperion's transport contract");

    // Central target (isocentre) vs an off-axis OAR, both inside the water cylinder.
    let mid_z = (NZ as f64 - 1.0) * SPACING / 2.0;
    let ptv = Dvh::from_volume_masked(
        &dose,
        spherical_mask(grid, Point3::new(CENTRE_MM, CENTRE_MM, mid_z), 8.0),
    );
    let oar = Dvh::from_volume_masked(
        &dose,
        spherical_mask(grid, Point3::new(CENTRE_MM - 16.0, CENTRE_MM, mid_z), 6.0),
    );
    assert!(
        ptv.count() > 0 && oar.count() > 0,
        "masks must select voxels"
    );

    // Rotational convergence: the central target is hotter than the off-axis OAR.
    assert!(
        ptv.mean() > oar.mean(),
        "PTV mean {} must exceed OAR mean {}",
        ptv.mean().into_base(),
        oar.mean().into_base()
    );
    let ptv_geud = ptv.generalized_eud(1.0).expect("valid PTV response");
    let oar_geud = oar.generalized_eud(1.0).expect("valid OAR response");
    assert!(ptv_geud > oar_geud, "PTV gEUD must exceed OAR gEUD");

    // Outcome models are well-formed probabilities. With TCD50 set below the PTV
    // gEUD the target controls (TCP > 0.5); with TD50 above the OAR gEUD the OAR is
    // spared (NTCP < 0.5) — both by ratio, independent of the absolute dose scale.
    assert!(
        *oar_geud.as_base() > 0.0,
        "OAR (in-water) must receive some dose"
    );
    let tcp = ptv
        .tcp_logistic(1.0, ptv_geud * 0.8, 2.0)
        .expect("valid tumour-control response");
    let ntcp = oar
        .ntcp_lkb(1.0, oar_geud * 2.0, 0.3)
        .expect("valid complication response");
    assert!(
        (0.0..=1.0).contains(&tcp) && tcp > 0.5,
        "PTV TCP {tcp} should exceed 0.5"
    );
    assert!(
        (0.0..=1.0).contains(&ntcp) && ntcp < 0.5,
        "OAR NTCP {ntcp} should be below 0.5"
    );
}

// Small helper: a positive normalization dose for gamma (max voxel dose).
trait PeakDose {
    fn mean_top_dose(&self) -> f64;
}
impl PeakDose for Volume<f64> {
    fn mean_top_dose(&self) -> f64 {
        let [nx, ny, nz] = self.grid().dims();
        let mut peak = 0.0_f64;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    peak = peak.max(self.get(i, j, k).unwrap());
                }
            }
        }
        peak
    }
}
