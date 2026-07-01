//! Typed error surface for the Helios stack.
//!
//! Libraries in Helios return explicit typed errors (never a stringly-typed
//! catch-all that collapses distinct failure modes). The enum is
//! `#[non_exhaustive]` so higher layers add their own failure modes as they are
//! built without a breaking change to downstream matches.

/// The canonical error type for the Helios foundation layer.
///
/// Additional variants are introduced by higher layers (geometry, physics, GPU,
/// I/O) as those layers are implemented; matches on this enum must include a
/// wildcard arm because it is `#[non_exhaustive]`.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum HeliosError {
    /// A domain value failed validation at a construction/parse boundary.
    ///
    /// Carries the offending field name, the rejected value, and the invariant
    /// that was violated so the message names both the value and the reason.
    #[error("invalid {field}: {value} — {reason}")]
    InvalidDomainValue {
        /// Name of the domain quantity that failed validation.
        field: &'static str,
        /// The rejected numeric value.
        value: f64,
        /// The invariant that `value` violated.
        reason: &'static str,
    },

    /// A medical-image I/O operation (e.g. DICOM parse/decode/attribute read)
    /// failed. Carries a human-readable description of the failing step; the
    /// distinct variant keeps I/O failures from collapsing into validation errors.
    #[error("DICOM I/O error: {reason}")]
    Dicom {
        /// Description of the parse/decode/attribute-lookup failure.
        reason: String,
    },
}

/// Convenience alias for fallible Helios operations.
pub type Result<T> = core::result::Result<T, HeliosError>;
