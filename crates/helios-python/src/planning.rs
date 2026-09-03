use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

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
#[expect(
    clippy::needless_pass_by_value,
    reason = "PyO3 extracts sequence arguments into owned Vec<f64>; a borrowed slice is not an extractable argument type"
)]
pub(crate) fn optimize_beam_weights(
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
