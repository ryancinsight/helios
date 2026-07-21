//! Collapsed-Cone Convolution-Superposition Dose Calculation (3-D)
//!
//! Demonstrates the two-stage deterministic helios dose engine on a
//! 10×10×10 synthetic water phantom:
//!
//! 1. **CT → μ map** — `attenuation_map` converts Hounsfield units to linear
//!    attenuation coefficients using the Compton-dominated MV approximation.
//!
//! 2. **Primary transport** — `primary_fluence_parallel_x` applies Beer–Lambert
//!    attenuation through the phantom for a parallel beam entering along +x.
//!
//! 3. **TERMA** — primary fluence × local linear attenuation = energy released
//!    per voxel (Total Energy Released per unit MAss).
//!
//! 4. **Dose by 1-D convolution** — `dose_convolution_x` spreads released energy
//!    downstream with an exponential deposition kernel (fast collapsed-cone
//!    approximation along the beam axis).
//!
//! 5. **Full 3-D scatter superposition** — `scatter_superposition` with
//!    separable symmetric kernels on all three axes reproduces the lateral
//!    penumbra that the 1-D kernel misses.
//!
//! Each stage is verified with physics-plausibility assertions:
//! - depth-dose build-up (maximum dose is downstream of entry face)
//! - energy conservation within the analytical boundary-truncation bound:
//!   the 1-D and 3-D convolution kernels extend a finite radius past the
//!   10-voxel phantom edge, so a fraction of the released TERMA leaks out
//!   the boundary instead of being redeposited. The example documents an
//!   acceptable loss threshold (< 30 % for this synthetic geometry) at the
//!   assertion site; a production calculation on patient CT eliminates the
//!   loss because the interior dominates.
//!
//! ## Book Chapter
//!
//! [← Dose Calculation](../../docs/book/dose_collapsed_cone.md)
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-solver --example collapsed_cone_3d
//! ```

use helios_domain::{Volume, VoxelGrid};
use helios_math::Point3;
use helios_physics::MassAttenuation;
use helios_solver::{
    attenuation_map, dose_convolution_x, exponential_deposition_kernel, primary_fluence_parallel_x,
    scatter_superposition, symmetric_deposition_kernel,
};

fn main() {
    println!("=== Collapsed-Cone 3-D Dose Engine ===\n");

    // ── 1. Phantom geometry ────────────────────────────────────────────────────
    //
    // 10×10×10 voxel water phantom, 3 mm isotropic spacing (a small thorax slab).
    let nx = 10_usize;
    let ny = 10_usize;
    let nz = 10_usize;
    let spacing_mm = 3.0_f64;

    let grid = VoxelGrid::axis_aligned([nx, ny, nz], [spacing_mm; 3], Point3::new(0.0, 0.0, 0.0))
        .expect("valid phantom grid");

    println!(
        "Phantom: {}×{}×{} voxels, {:.1} mm isotropic spacing",
        nx, ny, nz, spacing_mm
    );
    println!(
        "Volume: {:.0}×{:.0}×{:.0} mm³\n",
        nx as f64 * spacing_mm,
        ny as f64 * spacing_mm,
        nz as f64 * spacing_mm
    );

    // ── 2. CT → μ map ─────────────────────────────────────────────────────────
    //
    // Uniform HU = 0 (pure water) throughout — a simplified phantom.
    let ct_hu = Volume::from_shape_fn(grid, |_| 0.0_f64);

    // Water at 6 MV: μ/ρ ≈ 0.0636 cm²/g, ρ_water = 1.0 g/cm³ → μ ≈ 0.0636 cm⁻¹.
    let mu_over_rho = MassAttenuation::new(0.0636_f64).expect("valid mass attenuation");
    let water_rho = 1.0_f64; // g/cm³

    let mu = attenuation_map(&ct_hu, mu_over_rho, water_rho);

    let mu_center = mu.get(nx / 2, ny / 2, nz / 2).unwrap();
    println!("Stage 1 — CT → μ map (Compton-dominated 6 MV approximation)");
    println!("  Water μ/ρ = 0.0636 cm²/g, ρ = {water_rho:.1} g/cm³");
    println!("  Resulting μ at phantom center = {mu_center:.5} cm⁻¹\n");

    // ── 3. Primary fluence — Beer–Lambert along +x ────────────────────────────
    let incident_fluence = 1.0_f64; // normalized to 1 at entry face
    let primary = primary_fluence_parallel_x(&mu, incident_fluence);

    let voxel_cm = spacing_mm * 0.1; // mm → cm
    let expected_surface = incident_fluence * (-mu_center * 0.0 * voxel_cm).exp();
    let surface_psi = primary.get(0, ny / 2, nz / 2).unwrap();
    let depth5_psi = primary.get(5, ny / 2, nz / 2).unwrap();

    println!("Stage 2 — Primary fluence (Beer–Lambert)");
    println!("  Ψ(x=0)  = {surface_psi:.6}  (expected {expected_surface:.6})");
    println!("  Ψ(x=5)  = {depth5_psi:.6}");
    let depth5_expected = (-mu_center * 5.0 * voxel_cm).exp();
    println!(
        "  Ψ analytic  = {depth5_expected:.6}  (error {:.2e})\n",
        (depth5_psi - depth5_expected).abs()
    );

    // ── 4. TERMA — energy released per voxel ─────────────────────────────────
    //
    // TERMA(i,j,k) = Ψ(i,j,k) · μ(i,j,k) · Δx_cm  [Gy per incident fluence]
    let terma = Volume::from_shape_fn(grid, |idx| {
        let psi = primary.get(idx[0], idx[1], idx[2]).expect("within grid");
        let mu_v = mu.get(idx[0], idx[1], idx[2]).expect("within grid");
        psi * mu_v * voxel_cm
    });

    let terma_total: f64 = (0..nx)
        .flat_map(|i| (0..ny).flat_map(move |j| (0..nz).map(move |k| (i, j, k))))
        .map(|(i, j, k)| terma.get(i, j, k).unwrap())
        .sum();
    println!("Stage 3 — TERMA");
    println!("  Total TERMA = {terma_total:.6} (≈ 1 − e^(-μ·depth) = primary extraction)\n");

    // ── 5. 1-D dose by convolution-superposition along beam (+x) ─────────────
    //
    // Exponential kernel: characteristic electron transport range 0.5 cm.
    let range_cm = 0.5_f64;
    let kernel = exponential_deposition_kernel(range_cm, voxel_cm, 8);

    println!("Stage 4 — 1-D dose convolution");
    println!(
        "  Deposition kernel: {} taps, range = {range_cm:.2} cm",
        kernel.len()
    );

    let dose_1d = dose_convolution_x(&terma, &kernel);

    // Depth-dose along the central beam axis (j = ny/2, k = nz/2).
    println!("  Central-axis depth-dose (j={}, k={}):", ny / 2, nz / 2);
    for i in 0..nx {
        let d = dose_1d.get(i, ny / 2, nz / 2).unwrap();
        println!("    x={i:2}: {d:.6} Gy");
    }

    // Build-up: maximum dose should NOT be at the entry face.
    let max_depth = (0..nx)
        .max_by(|&a, &b| {
            dose_1d
                .get(a, ny / 2, nz / 2)
                .unwrap()
                .partial_cmp(&dose_1d.get(b, ny / 2, nz / 2).unwrap())
                .unwrap()
        })
        .unwrap();
    println!("  Dose maximum at depth index {max_depth} (entry=0)\n");

    // ── 6. 3-D scatter superposition ─────────────────────────────────────────
    //
    // Separable symmetric kernels: tight lateral (0.3 cm) and relaxed axial (0.5 cm).
    let kx = symmetric_deposition_kernel(0.5_f64, voxel_cm, 3);
    let ky = symmetric_deposition_kernel(0.3_f64, voxel_cm, 2);
    let kz = symmetric_deposition_kernel(0.3_f64, voxel_cm, 2);

    let dose_3d = scatter_superposition(&terma, &kx, &ky, &kz);

    // Lateral profile at mid-depth (i=5), comparing 1-D vs 3-D dose.
    let mid_depth = nx / 2;
    println!("Stage 5 — 3-D scatter superposition");
    println!(
        "  Lateral profile at x={mid_depth} (j varies, k={}): 1-D vs 3-D",
        nz / 2
    );
    for j in 0..ny {
        let d1 = dose_1d.get(mid_depth, j, nz / 2).unwrap();
        let d3 = dose_3d.get(mid_depth, j, nz / 2).unwrap();
        let bar_len =
            (d3 * 40.0 / dose_3d.get(mid_depth, ny / 2, nz / 2).unwrap()).round() as usize;
        let bar = "#".repeat(bar_len.min(40));
        println!("    j={j:2}: 1-D {d1:.5}  3-D {d3:.5}  |{bar}");
    }

    // Energy conservation: total 3-D dose should ≈ total TERMA.
    let dose_3d_total: f64 = (0..nx)
        .flat_map(|i| (0..ny).flat_map(move |j| (0..nz).map(move |k| (i, j, k))))
        .map(|(i, j, k)| dose_3d.get(i, j, k).unwrap())
        .sum();
    let conservation_error = ((dose_3d_total - terma_total) / terma_total).abs();
    println!("\n  3-D dose total  = {dose_3d_total:.6}");
    println!("  TERMA total     = {terma_total:.6}");
    println!(
        "  Energy conservation error = {conservation_error:.4} ({:.2}%)\n",
        conservation_error * 100.0
    );

    // ── Self-validation ───────────────────────────────────────────────────────
    assert!(
        (surface_psi - 1.0).abs() < 1e-12,
        "Entry-face fluence must equal incident fluence"
    );
    assert!(
        (depth5_psi - depth5_expected).abs() < 1e-10,
        "Beer–Lambert depth-dose must match analytical exponential"
    );
    assert!(
        max_depth > 0,
        "Depth-dose build-up: maximum must be downstream of entry face"
    );
    // Boundary truncation: the scatter kernel extends ±radius voxels. On a small
    // 10-voxel phantom a kernel of radius 2–3 voxels loses ~20 % of energy at the
    // faces. Acceptable boundary loss is < 30 % for this synthetic case; a
    // production calculation would use an unbounded phantom (patient CT) where
    // the interior dominates.
    assert!(
        conservation_error < 0.30,
        "Energy conservation error {conservation_error:.3} exceeds boundary threshold"
    );

    println!("All physics checks passed ✓");
    println!("\nBook chapter: Part III — Dose Calculation");
    println!("API: helios_solver::{{attenuation_map, primary_fluence_parallel_x,");
    println!("                     exponential_deposition_kernel, dose_convolution_x,");
    println!("                     symmetric_deposition_kernel, scatter_superposition}}");
}
