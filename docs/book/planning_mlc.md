# Chapter 13 — MLC Models and Leaf Sequencing

The Multi-Leaf Collimator (MLC) is the binary modulator in a TomoTherapy
machine. Each of the 64 binary leaves (1.0 cm width at isocentre) is either
open or closed during delivery.

`
ust
use helios_domain::{LeafOpenTimeSinogram, MlcModel};

let mlc = MlcModel::binary(64, 1.0); // 64 leaves, 1 cm width
let lots = LeafOpenTimeSinogram { leaves: 64, gantry_angles: 51, data: ... };
`

## Leaf Open Time

LeafOpenTimeSinogram (LOTS) encodes the fraction of each gantry rotation
that each leaf is open: values in [0, 1]. At 0 the leaf is always closed,
at 1 always open.

## Tongue-and-Groove Effect

Adjacent leaf pairs share a slight overlap (tongue-and-groove feature) that
reduces interleaf leakage. MlcModel accounts for this in fluence calculation.

## Further Reading

- [Helical Delivery and Sinogram](planning_helical.md)
- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)