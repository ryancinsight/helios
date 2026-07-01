//! Helios time-dependent helical delivery / acquisition simulation.
//!
//! Combines the [`HelicalDelivery`](helios_domain::HelicalDelivery) kinematics
//! (gantry rotation + couch translation) with the [`helios_solver`] forward
//! projector to simulate a helical scan over time: at each projection the gantry
//! angle rotates the beam in the axial plane while the couch advances the imaged
//! slice along the rotation axis, tracing a helix through the patient. Each
//! projection's central ray is forward-projected through the attenuation volume to
//! yield the MVCT detector signal (optical depth and transmission).
//!
//! This is the CPU reference; parallel/GPU orchestration (moirai/hephaestus over
//! independent projections) is a later increment. Grids are axis-aligned
//! (`helios-domain` `VoxelGrid`).
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod acquisition;
mod delivery;
mod dose_accumulation;

pub use acquisition::{simulate_helical_sinogram, HelicalProjection};
pub use delivery::{simulate_helical_delivery, total_delivered_fluence, DeliveryFrame};
pub use dose_accumulation::accumulate_delivered_dose;
