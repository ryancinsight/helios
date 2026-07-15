# Parallel-Beam Radon Transform

The Radon transform forward-projects a 2-D attenuation map μ(x, y) along
parallel rays at 
_angles uniformly-spaced rotation angles:

`	ext
R[θ, t] = ∫ μ(t·cos θ − s·sin θ, t·sin θ + s·cos θ) ds
`

In Helios the Radon transform is provided by helios-imaging:

`ust
use helios_imaging::parallel_beam_radon;

let n_angles = 180;
let sinogram = parallel_beam_radon(&mu_volume, n_angles);
// sinogram.grid().dims() == [n_angles, nx]
`

## Physical Interpretation

Each row of the sinogram is the projection at one gantry angle.
For a cylinder of uniform attenuation μ₀:

`	ext
R[θ, 0] = μ₀ · 2√(r² − t²)  for |t| < r
`

This analytical oracle is used in the regression tests.

## Performance

The implementation is parallelised via moirai::map_collect_index_with —
each angle is computed independently on a separate thread.

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Filtered Back Projection](imaging_fbp.md)
- [Dose Calculation: Terma](dose_terma.md)
