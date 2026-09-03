//! Thin `PyO3` binding surface for Helios (`import helios`).
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
#![deny(missing_docs)]

mod physics;
mod planning;

use physics::{
    compton_mass_attenuation, klein_nishina_cross_section, mass_density_from_hu,
    thomson_cross_section,
};
use planning::optimize_beam_weights;

use pyo3::prelude::*;

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
