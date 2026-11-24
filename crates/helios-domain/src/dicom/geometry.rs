//! DICOM orientation invariants and the local rotation basis.

use helios_math::{Scalar, UnitQuaternion, Vector3};
use helios_quantities::HeliosError;

pub(super) const ORIENTATION_TOL: f64 = 1.0e-6;

pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(super) fn normalize3(v: [f64; 3], component: &'static str) -> Result<[f64; 3], HeliosError> {
    let norm2 = dot3(v, v);
    if !norm2.is_finite() || norm2 <= ORIENTATION_TOL * ORIENTATION_TOL {
        return Err(HeliosError::Dicom {
            reason: format!("{component} is zero-length or non-finite"),
        });
    }
    let inv = norm2.sqrt().recip();
    Ok([v[0] * inv, v[1] * inv, v[2] * inv])
}

pub(super) fn rotation_from_axes<T: Scalar>(
    row_dir: [f64; 3],
    col_dir: [f64; 3],
    normal_dir: [f64; 3],
) -> Result<UnitQuaternion<T>, HeliosError> {
    UnitQuaternion::try_from_rotation_columns(
        Vector3::new(
            T::from_f64(row_dir[0]),
            T::from_f64(row_dir[1]),
            T::from_f64(row_dir[2]),
        ),
        Vector3::new(
            T::from_f64(col_dir[0]),
            T::from_f64(col_dir[1]),
            T::from_f64(col_dir[2]),
        ),
        Vector3::new(
            T::from_f64(normal_dir[0]),
            T::from_f64(normal_dir[1]),
            T::from_f64(normal_dir[2]),
        ),
        T::from_f64(ORIENTATION_TOL),
    )
    .map_err(|e| HeliosError::Dicom {
        reason: format!("invalid ImageOrientationPatient basis: {e}"),
    })
}

#[inline]
pub(super) fn allclose3(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    (0..3).all(|axis| (a[axis] - b[axis]).abs() <= tol)
}
