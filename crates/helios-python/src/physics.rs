use aequitas::systems::si::{
    quantities::Energy,
    units::{MegaElectronVolt, SquareCentimeterPerGram},
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Reject a photon energy that is not strictly positive and finite (`MeV`).
pub(crate) fn checked_energy_mev(energy_mev: f64) -> PyResult<f64> {
    if !energy_mev.is_finite() || energy_mev <= 0.0 {
        return Err(PyValueError::new_err(
            "photon energy must be a finite, strictly positive value in MeV",
        ));
    }
    Ok(energy_mev)
}

/// Total Thomson (classical) scattering cross-section `σ_T` (m²/electron).
#[pyfunction]
pub(crate) fn thomson_cross_section() -> f64 {
    helios_physics::thomson_cross_section::<f64>()
}

/// Total Klein–Nishina Compton cross-section (m²/electron) at `energy_mev`.
///
/// Raises `ValueError` if `energy_mev` is not finite and positive.
#[pyfunction]
pub(crate) fn klein_nishina_cross_section(energy_mev: f64) -> PyResult<f64> {
    let e = checked_energy_mev(energy_mev)?;
    Ok(helios_physics::klein_nishina_cross_section::<f64>(
        Energy::from_unit::<MegaElectronVolt>(e),
    ))
}

/// Compton mass attenuation coefficient (μ/ρ, cm²/g) at `energy_mev` for a
/// material of effective ⟨Z/A⟩ `z_over_a` (water ≈ 0.5551), derived as
/// `σ_KN(E) · (electrons per gram)`.
///
/// Raises `ValueError` if `energy_mev` is not finite and positive.
#[pyfunction]
pub(crate) fn compton_mass_attenuation(energy_mev: f64, z_over_a: f64) -> PyResult<f64> {
    let e = checked_energy_mev(energy_mev)?;
    let electrons_per_gram = helios_physics::electrons_per_gram::<f64>(z_over_a);
    Ok(helios_physics::compton_mass_attenuation::<f64>(
        Energy::from_unit::<MegaElectronVolt>(e),
        electrons_per_gram,
    )
    .in_unit::<SquareCentimeterPerGram>())
}

/// Mass density (g/cm³) from a Hounsfield unit via bilinear CT calibration,
/// given the reference `water_density_g_cm3`.
#[pyfunction]
pub(crate) fn mass_density_from_hu(hu: f64, water_density_g_cm3: f64) -> f64 {
    helios_physics::mass_density_from_hu::<f64>(hu, water_density_g_cm3)
}
