//! End-to-end helical TomoTherapy + MVCT workflow demo (H-041b).
//!
//! Runs the full Helios pipeline on a synthetic CT phantom and renders inspectable
//! PNG artifacts (per the Output & visual verification discipline):
//!
//! 1. CT phantom (HU) → attenuation map `μ`.
//! 2. Imaging: parallel-beam Radon → FBP reconstruction of `μ`.
//! 3. Therapy: helical MLC delivery → divergent-fan terma → beam-following
//!    poly-energetic (beam-hardened) collapsed-cone dose; DVH + 3%/2 mm self-gamma.
//!
//! Writes `ct.png`, `mu.png`, `recon.png`, `dose.png` to the output directory (first
//! CLI arg, default `helios_workflow_output/`) and prints a quantitative summary.
//!
//! Run: `cargo run -p helios-simulation --example tomotherapy_workflow [out_dir]`.

use std::path::{Path, PathBuf};

use aequitas::systems::si::{
    quantities::AreaPerMass,
    units::{Gray, SquareCentimeterPerGram},
};
use helios_analysis::{gamma_index_3d, gamma_pass_rate, roi_statistics, Dvh};
use helios_domain::{HelicalDelivery, LeafOpenTimeSinogram, MlcModel, Volume, VoxelGrid};
use helios_imaging::{filtered_back_projection, parallel_beam_radon};
use helios_math::Point3;
use helios_simulation::{
    accumulate_delivered_dose_anisotropic, simulate_helical_delivery, BeamGeometry, CollapsedCone,
    SpectralComponent,
};
use helios_solver::attenuation_map;
use hyperion::coefficient::MassAttenuation;

const NX: usize = 61;
const NZ: usize = 5;
const SPACING: f64 = 1.0;

fn centre_mm() -> f64 {
    (NX as f64 - 1.0) * SPACING / 2.0
}

/// Air / water-cylinder / bone-insert CT phantom (HU).
fn ct_phantom() -> Volume<f64> {
    let grid =
        VoxelGrid::axis_aligned([NX, NX, NZ], [SPACING; 3], Point3::new(0.0, 0.0, 0.0)).unwrap();
    let c = centre_mm();
    Volume::from_shape_fn(grid, move |idx| {
        let dx = idx[0] as f64 * SPACING - c;
        let dy = idx[1] as f64 * SPACING - c;
        let r = (dx * dx + dy * dy).sqrt();
        if r <= 8.0 {
            800.0
        } else if r <= 25.0 {
            0.0
        } else {
            -1000.0
        }
    })
}

/// Render slice `k` of `vol` to a normalized grayscale PNG (max → white).
fn render_slice(vol: &Volume<f64>, k: usize, path: &Path) {
    let [nx, ny, _] = vol.grid().dims();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for j in 0..ny {
        for i in 0..nx {
            let v = vol.get(i, j, k).unwrap();
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    let span = if hi > lo { hi - lo } else { 1.0 };
    let mut img = image::GrayImage::new(nx as u32, ny as u32);
    for j in 0..ny {
        for i in 0..nx {
            let v = vol.get(i, j, k).unwrap();
            let g = (((v - lo) / span) * 255.0).round().clamp(0.0, 255.0) as u8;
            img.put_pixel(i as u32, j as u32, image::Luma([g]));
        }
    }
    img.save(path).expect("write PNG");
}

fn peak(vol: &Volume<f64>) -> f64 {
    let [nx, ny, nz] = vol.grid().dims();
    let mut p = 0.0_f64;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                p = p.max(vol.get(i, j, k).unwrap());
            }
        }
    }
    p
}

fn main() {
    let out: PathBuf = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "helios_workflow_output".to_owned())
        .into();
    std::fs::create_dir_all(&out).expect("create output dir");

    // 1. CT → μ.
    let ct = ct_phantom();
    let coefficient = MassAttenuation::new(AreaPerMass::from_unit::<SquareCentimeterPerGram>(0.06))
        .expect("valid water mass attenuation");
    let mu = attenuation_map(&ct, coefficient, 1.0).expect("fixture calibration is finite");

    // 2. Imaging: Radon → FBP.
    let angles: Vec<f64> = (0..180)
        .map(|a| a as f64 * std::f64::consts::PI / 180.0)
        .collect();
    let n_off = 121;
    let offsets: Vec<f64> = (0..n_off)
        .map(|j| -35.0 + j as f64 * 70.0 / (n_off - 1) as f64)
        .collect();
    let sino = parallel_beam_radon(&mu, &angles, &offsets, 500.0, 0.25);
    let recon_grid =
        VoxelGrid::axis_aligned([NX, NX, 1], [SPACING; 3], Point3::new(0.0, 0.0, 0.0)).unwrap();
    let recon = filtered_back_projection(&sino, &recon_grid);
    let water = roi_statistics(&recon, [40, 28, 0], [45, 33, 1]);

    // 3. Therapy: helical MLC delivery → dose.
    let (leaves, projections) = (21, 40);
    let delivery = HelicalDelivery::new(51, 20.0, 0.3, 10.0, 0.0, 0.0).unwrap();
    // Modulated fluence: a central band of open leaves (a simple target aperture).
    let mut fractions = vec![0.0; projections * leaves];
    for p in 0..projections {
        for l in 6..15 {
            fractions[p * leaves + l] = 1.0;
        }
    }
    let lot = LeafOpenTimeSinogram::from_fractions(projections, leaves, fractions).unwrap();
    let mlc = MlcModel::new(0.01, 0.05).unwrap();
    let frames = simulate_helical_delivery(&delivery, &lot, &mlc);
    // Beam-following, poly-energetic collapsed cone: a two-component (soft/hard)
    // spectrum, forward-peaked (longer downstream range) and re-oriented to each
    // frame's gantry direction. Beam hardening + forward transport shape the dose.
    let voxel_cm = SPACING / 10.0;
    let spectrum = [
        SpectralComponent {
            range_up_cm: 0.10,
            range_down_cm: 0.35,
            weight: 0.7,
        }, // soft component
        SpectralComponent {
            range_up_cm: 0.05,
            range_down_cm: 0.90,
            weight: 0.3,
        }, // hard component (reaches farther downstream)
    ];
    let cone = CollapsedCone::poly_forward_peaked(&spectrum, 0.6, voxel_cm, 2, 3, 2);
    let dose = accumulate_delivered_dose_anisotropic(
        &frames,
        &mu,
        BeamGeometry::PointSource {
            source_axis_mm: 850.0,
        },
        3.0,
        0.25,
        &cone,
    )
    .expect("attenuation map satisfies Hyperion's transport contract");

    // Analysis.
    let dvh = Dvh::from_volume(&dose);
    let peak_dose = peak(&dose);
    let gamma = gamma_index_3d(&dose, &dose, 0.03, 2.0, peak_dose, 6.0).unwrap();
    let pass = gamma_pass_rate(&gamma, &dose, 0.1 * peak_dose);

    // Renders.
    render_slice(&ct, NZ / 2, &out.join("ct.png"));
    render_slice(&mu, NZ / 2, &out.join("mu.png"));
    render_slice(&recon, 0, &out.join("recon.png"));
    render_slice(&dose, NZ / 2, &out.join("dose.png"));

    println!("Helios end-to-end workflow ({}³ voxels)", NX);
    println!(
        "  μ:      water 0.0600 cm⁻¹, bone {:.4}, air ~0",
        mu.get(NX / 2, NX / 2, NZ / 2).unwrap()
    );
    println!(
        "  recon:  water ROI μ = {:.4} cm⁻¹ (target 0.0600, {:+.1}%)",
        water.mean,
        (water.mean / 0.06 - 1.0) * 100.0
    );
    println!(
        "  dose:   total {:.3}, peak {:.4}, DVH mean {:.4}",
        dose.sum(),
        peak_dose,
        dvh.mean().in_unit::<Gray>()
    );
    println!(
        "  gamma:  3%/2 mm self-consistency pass rate {:.1}%",
        pass * 100.0
    );
    println!(
        "  wrote ct.png, mu.png, recon.png, dose.png to {}",
        out.display()
    );
}
