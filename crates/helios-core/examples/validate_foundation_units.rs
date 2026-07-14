//! Helios-core example: typestate slots for the foundation layer.
//!
//! Exercises the validating newtypes `EnergyMeV`, `HounsfieldUnit`,
//! `VoxelSpacingMm` — the typed slots every higher layer passes at construction
//! boundaries. Establishes that the typestate slots reject everything the physics
//! layer must not see (non-finite values, out-of-range CT numbers, non-positive
//! spacings), so upstream callers must build domain values from validated sources
//! rather than hand-poking the underlying float.
//!
//! Run with:  cargo run --example validate_foundation_units -p helios-core

use helios_core::{EnergyMeV, HeliosError, HounsfieldUnit, VoxelSpacingMm};

/// Returns the canonical clinical water slot: 6 MV beam, 0 HU water, 1.0 mm
/// voxel pitch. Demonstrates that the three slots carry independent unit
/// domains and round-trip identity per axis.
fn canonical_water_slot() -> (EnergyMeV, HounsfieldUnit, VoxelSpacingMm) {
    let energy = EnergyMeV::try_from(6.0).expect("6 MeV is a clinically valid beam energy");
    let hu = HounsfieldUnit::try_from(0.0).expect("0 HU is the water CT number by convention");
    let spacing = VoxelSpacingMm::try_from(1.0).expect("1.0 mm voxel spacing is strictly positive");
    (energy, hu, spacing)
}

/// Asserts the typestate reject paths — every value that should not survive
/// boundary validation must surface a `HeliosError::InvalidDomainValue`.
fn assert_typestate_rejects() {
    fn invalid<T>(result: Result<T, HeliosError>) {
        match result {
            Err(HeliosError::InvalidDomainValue { .. }) => {}
            Err(other) => panic!("expected InvalidDomainValue, got {other:?}"),
            Ok(_) => panic!("expected Err, got Ok"),
        }
    }

    // Energy: non-finite, zero, and negative must each be rejected.
    invalid(EnergyMeV::try_from(0.0));
    invalid(EnergyMeV::try_from(-2.5));
    invalid(EnergyMeV::try_from(f64::NAN));
    invalid(EnergyMeV::try_from(f64::INFINITY));

    // Hounsfield: range is calibrated to the CT-number scale [-1024, 31743].
    invalid(HounsfieldUnit::try_from(-1025.0));
    invalid(HounsfieldUnit::try_from(HounsfieldUnit::MAX + 1.0));

    // Spacing: positivity is the only carrier-invariant.
    invalid(VoxelSpacingMm::try_from(0.0));
    invalid(VoxelSpacingMm::try_from(-0.5));

    // Boundary values must succeed (inclusive ranges).
    EnergyMeV::try_from(f64::MIN_POSITIVE).expect("MIN_POSITIVE succeeds");
    HounsfieldUnit::try_from(HounsfieldUnit::MIN).expect("MIN is inclusive");
    HounsfieldUnit::try_from(HounsfieldUnit::MAX).expect("MAX is inclusive");
    VoxelSpacingMm::try_from(f64::MIN_POSITIVE).expect("MIN_POSITIVE succeeds");
}

fn main() {
    let (energy, hu, spacing) = canonical_water_slot();
    println!("Canonical water slot: energy={energy}, hu={hu}, spacing={spacing}");

    // Display carries the unit suffix by construction.
    assert_eq!(energy.to_string(), "6 MeV");
    assert_eq!(hu.to_string(), "0 HU");
    assert_eq!(spacing.to_string(), "1 mm");

    // Exact round-trip identity by axis.
    assert_eq!(energy.get(), 6.0);
    assert_eq!(hu.get(), 0.0);
    assert_eq!(spacing.get(), 1.0);

    assert_typestate_rejects();
    println!("All typestate slots round-trip and reject every invalid input.");
}
