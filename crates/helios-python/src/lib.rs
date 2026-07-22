//! Thin PyO3 binding surface for Helios (`import helios`).
//!
//! This is the ONLY Helios crate permitted to depend on `pyo3`. It holds no
//! domain logic: every function validates and converts its Python arguments into
//! typed Rust values, calls the corresponding `helios-*` core function, maps
//! [`helios_physics`]/[`helios_planning`] failures to a Python `ValueError`, and
//! converts the result back. Compute-heavy calls release the GIL via
//! [`Python::detach`] so Python threads run concurrently with the Rust
//! core. Concrete `f64` is used at this FFI boundary (the sanctioned place for a
//! concrete numeric type); the underlying kernels remain generic over `Scalar`.
#![forbid(unsafe_code)]

use aequitas::systems::si::units::SquareCentimeterPerGram;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Reject a photon energy that is not strictly positive and finite (MeV).
fn checked_energy_mev(energy_mev: f64) -> PyResult<f64> {
    if !energy_mev.is_finite() || energy_mev <= 0.0 {
        return Err(PyValueError::new_err(
            "photon energy must be a finite, strictly positive value in MeV",
        ));
    }
    Ok(energy_mev)
}

/// Total Thomson (classical) scattering cross-section σ_T (m²/electron).
#[pyfunction]
fn thomson_cross_section() -> f64 {
    helios_physics::thomson_cross_section::<f64>()
}

/// Total Klein–Nishina Compton cross-section (m²/electron) at `energy_mev`.
///
/// Raises `ValueError` if `energy_mev` is not finite and positive.
#[pyfunction]
fn klein_nishina_cross_section(energy_mev: f64) -> PyResult<f64> {
    let e = checked_energy_mev(energy_mev)?;
    Ok(helios_physics::klein_nishina_cross_section::<f64>(e))
}

/// Compton mass attenuation coefficient (μ/ρ, cm²/g) at `energy_mev` for a
/// material of effective ⟨Z/A⟩ `z_over_a` (water ≈ 0.5551), derived as
/// `σ_KN(E) · (electrons per gram)`.
///
/// Raises `ValueError` if `energy_mev` is not finite and positive.
#[pyfunction]
fn compton_mass_attenuation(energy_mev: f64, z_over_a: f64) -> PyResult<f64> {
    let e = checked_energy_mev(energy_mev)?;
    let electrons_per_gram = helios_physics::electrons_per_gram::<f64>(z_over_a);
    Ok(
        helios_physics::compton_mass_attenuation::<f64>(e, electrons_per_gram)
            .in_unit::<SquareCentimeterPerGram>(),
    )
}

/// Mass density (g/cm³) from a Hounsfield unit via bilinear CT calibration,
/// given the reference `water_density_g_cm3`.
#[pyfunction]
fn mass_density_from_hu(hu: f64, water_density_g_cm3: f64) -> f64 {
    helios_physics::mass_density_from_hu::<f64>(hu, water_density_g_cm3)
}

/// Projected-gradient inverse-planning optimum: non-negative beam weights
/// minimizing `½‖A·x − d‖²`.
///
/// `influence` is the row-major dose-influence matrix `A` of shape
/// `voxels × beamlets`; `prescription` is the target dose `d` (length `voxels`).
/// The GIL is released around the iterative solve.
///
/// Raises `ValueError` if `influence` length ≠ `voxels·beamlets` or
/// `prescription` length ≠ `voxels`.
#[pyfunction]
#[pyo3(signature = (influence, voxels, beamlets, prescription, iterations, step))]
fn optimize_beam_weights(
    py: Python<'_>,
    influence: Vec<f64>,
    voxels: usize,
    beamlets: usize,
    prescription: Vec<f64>,
    iterations: usize,
    step: f64,
) -> PyResult<Vec<f64>> {
    let dose_influence = helios_planning::DoseInfluence::from_rows(voxels, beamlets, influence)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    if prescription.len() != voxels {
        return Err(PyValueError::new_err(format!(
            "prescription length {} does not match voxel count {voxels}",
            prescription.len()
        )));
    }
    let weights = py.detach(|| {
        helios_planning::optimize_beam_weights(&dose_influence, &prescription, iterations, step)
    });
    Ok(weights)
}

/// Helios Python module.
#[pymodule]
fn helios(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(thomson_cross_section, m)?)?;
    m.add_function(wrap_pyfunction!(klein_nishina_cross_section, m)?)?;
    m.add_function(wrap_pyfunction!(compton_mass_attenuation, m)?)?;
    m.add_function(wrap_pyfunction!(mass_density_from_hu, m)?)?;
    m.add_function(wrap_pyfunction!(optimize_beam_weights, m)?)?;
    Ok(())
}
