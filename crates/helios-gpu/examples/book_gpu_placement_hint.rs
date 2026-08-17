//! Book example: GPU attenuation mapping through provider-owned placement hints.
//!
//! `GpuAttenuationMapper` routes uploads/allocations through Themis placement
//! hints internally (`PlacementHint::Tier(MemoryTier::Device)`), so callers keep
//! a simple physics-facing API.

use aequitas::systems::si::{
    quantities::{AreaPerMass, MassDensity},
    units::{GramPerCubicCentimeter, SquareCentimeterPerGram},
};
use helios_gpu::{default_device, GpuAttenuationMapper};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let device = match default_device() {
        Ok(device) => device,
        Err(_) => return Ok(()),
    };

    let mapper = GpuAttenuationMapper::new(
        device,
        AreaPerMass::from_unit::<SquareCentimeterPerGram>(0.0636_f32),
        MassDensity::from_unit::<GramPerCubicCentimeter>(1.0_f32),
    )?;

    let hu = vec![0.0_f32, 150.0, 700.0, -100.0];
    let mut mu = vec![0.0_f32; hu.len()];
    mapper.map_into(&hu, &mut mu)?;

    let expected = hu
        .iter()
        .map(|v| (0.0636_f32 * (1.0 + v / helios_core::constants::HU_SCALE_DENOMINATOR as f32)).max(0.0))
        .collect::<Vec<_>>();

    for (gpu, cpu_like) in mu.iter().zip(expected.iter()) {
        assert!(
            (gpu - cpu_like).abs() < 1e-4,
            "gpu={gpu} expected={cpu_like}"
        );
    }

    Ok(())
}
