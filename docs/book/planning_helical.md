# Chapter 14 — Helical Delivery and Sinogram

TomoTherapy delivers dose helically: the gantry rotates while the couch
translates, producing a helical trajectory of the beam relative to the patient.

`
ust
use helios_domain::{HelicalDelivery, LeafOpenTimeSinogram};

let delivery = HelicalDelivery {
    sinogram: lots,
    beam_geometry: BeamGeometry {
        gantry_period_s: 20.0,   // seconds per rotation
        pitch: 0.287,             // couch advance per rotation / field width
        field_width_cm: 2.5,
    },
};
`

## Pitch

Pitch = (couch advance per rotation) / (field width). Typical values: 0.2–0.4.
Lower pitch = more overlap = higher dose uniformity but longer treatment.

## Sinogram

The delivery sinogram has shape [n_gantry_angles, n_leaves].
For a standard TomoTherapy treatment: 51 gantry projections × 64 leaves.

## Further Reading

- [MLC Models and Leaf Sequencing](planning_mlc.md)
- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
