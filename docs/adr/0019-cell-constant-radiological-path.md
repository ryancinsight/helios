# ADR 0019: Cell-constant radiological path beside the trilinear march

- Status: Proposed
- Date: 2026-11-13
- Scope: `helios-solver` projection; board item H-011d

## Context

`helios_solver::forward_project_ray(mu, ray, step_mm)` integrates the
node-centred trilinear field. Integer `VoxelGrid` indices are voxel-centre
nodes. The ray is clipped to the node box `[0, n-1]` in index space and the
chord is midpoint-sampled with `Volume::sample_trilinear`. `GpuProjector` and
the Hephaestus `RayIntegralOps` kernels (locked revision `72ffa05`) integrate
the same model and are differentially tested against it. TERMA deposition
shares `RayMarchPlan`.

H-011d asks for an exact path, which exists only for a field model whose line
integral has a closed form. According to the abstract retrieved from PubMed
(Siddon 1985, Med. Phys. 12(2):252-255, PMID 4000088,
[DOI 10.1118/1.595715](https://doi.org/10.1118/1.595715)), the exact method
treats the CT array as the intersection volumes of three orthogonal sets of
equally spaced parallel planes, so its cost scales with the 3N planes rather
than the N³ voxels. The paper body was not read; the formulas below are
derived here.

## Decision

### 1. Field model, support, interval

Each stored value is a constant attenuation over its cell: voxel `i` spans
`[i - 1/2, i + 1/2]` on each local index axis, and the support is the cell box
`[-1/2, n - 1/2]`. The values remain voxel-centre samples, so the grid
convention does not change. The support does change: the cell box extends half
a voxel past the node box on every face. For uniform `mu` the two projectors
therefore differ exactly by `mu (L_cell - L_node)`; the 11-node, 2 mm projector
fixture gives 22 mm against 20 mm. A one-node axis has a zero-width node box
but a full cell.

Tie rule: cells are half-open `[i - 1/2, i + 1/2)`, and the last cell on each
axis is closed. Every point of the closed box therefore belongs to exactly one
cell. A ray in an interior face belongs to the upper cell. A ray in the upper
outer face belongs to cell `n - 1`.

Interval: the half-line part `[t0, t1] = [max(t_enter, 0), t_exit]` of the ray
clipped to the cell box. This matches the Hephaestus contract at `72ffa05`
(host `volume.rs:53`, wgpu `volume.rs:99`) and the gaia `Ray` half-line. A miss
returns `None`, and a zero-length contact returns `Some(0)`, as the march does.
The CPU march and TERMA deposition integrate the same half-line: the shared
grid interval clamps the signed slab entry to the origin.

### 2. Traversal

Map the ray into the grid's local millimetre frame with the inverse pose. The
pose is an isometry, so `d` stays unit and `t` stays in world mm. Let
`p = o + t0 d` be the entry point. Plane `k` of axis `a` lies at
`x = (k - 1/2) s_a`. Each crossing is measured from the entry point:

```text
r_a(k) = ((k - 1/2) s_a - p_a) / d_a    in [0, L], L = t1 - t0, d_a != 0
```

Each `r` is computed from its plane index, so rounding does not accumulate along
the ray. An axis with `d_a == 0` has no crossings, and the tie rule fixes its
index from `p_a`. Merging the three sorted crossing sequences with `{0, L}` and
removing equal values gives the segment boundaries, so a vertex crossing yields
no zero-length segment. Each segment lies in the cell containing its midpoint
and has length `r_{j+1} - r_j`. There are at most `N <= nx + ny + nz - 2`
segments.

### 3. API

The API is additive, in a new `helios-solver` module beside the clip whose pose
transform it shares:

- `cell_path(grid, ray) -> Option<CellPath<T>>`: an allocation-free
  `Iterator<Item = CellSegment<T>> + FusedIterator`. Each segment carries
  `index: [usize; 3]` and `length: Length<T>` (Aequitas, per ADR 0016).
- `radiological_path(mu, ray) -> Option<T>`: `tau = sum(mu_k l_k)`, with `mu`
  in cm⁻¹, summed in traversal order and converted mm→cm once.

`forward_project_ray`, `RayMarchPlan`, and the GPU projector are unchanged.
Identifiers name the operation, not the algorithm's author.

## Alternatives rejected

- One projector parameterised by a field-model policy (the strongest
  alternative): `forward_project_ray::<M: FieldModel>` with
  `NodeTrilinear { step }` and `CellConstant`. Rejected because the two models
  share only the pose transform, which the clip helper already shares. They
  differ in support, in input (a step exists only for the march), and in output
  (per-cell lengths, which a matched projector/back-projector needs, cannot come
  from the march). Overturned by a call site that must select the model at run
  time; that seam then belongs on the consumer, for example a projector
  strategy taken by `parallel_beam_radon`.
- Replace the march: rejected. Every imaging, acquisition, portal, and dose
  result would shift by the border shell and the cell steps, and CPU/GPU parity
  would break until Hephaestus ships a cell-traversal kernel.
- Host the traversal in Gaia: rejected for now. A search of the Hephaestus and
  Gaia default branches for siddon, amanatides, and voxel/grid traversal found
  no second consumer; one promotes it.

## Error bound

Let `u = T::EPSILON / 2` and `γ_n = n u / (1 - n u)`. Define `E` as the largest
box extent `n_a s_a` (mm), `R0 = |o|_inf + |pose translation|_inf + |t0|`, and
`d_min` as the smallest nonzero `|d_a|`. On an identity pose each `r` carries
seven roundings proportional to `R0 + E`: three in the clip's `t0` (subtract,
reciprocal, multiply), one in the origin translation, two in `p`, and one in the
plane coordinate (`k - 1/2` is exact for `k` below `2^23` in `f32` and `2^52` in
`f64`). The numerator subtraction adds `u (E + R0)`, and the division adds
`u L`. With both clip endpoints treated as crossings, every boundary therefore
satisfies `|dr| <= eta = γ_7 (R0 + E) / d_min + u L`. An oriented pose adds
rotation rounding to `o` and `d`, which gives `eta_rot = γ_9 (R0 + E + L) / d_min`.

Each segment length then errs by at most `2 eta + u l`. The products, the sum of
`N` terms, and the mm→cm scaling (one rounding for `1/10`, one for the product)
contribute the rest. To first order:

```text
|tau_hat - tau| <= (mu_max / 10) (2 N eta + γ_{N+1} L) + γ_2 |tau|
```

Tests evaluate this bound in a named helper; there are no literal tolerances.
Test rays keep at least `eta` from cell edges, or lie exactly in a plane under
the tie rule, so the cell sequence is stable.

## Acceptance oracle (f32 and f64)

1. Uniform box: an axis ray gives `tau = mu n s / 10`. An oblique ray's segment
   lengths sum to the slab-clip chord of the cell box.
2. Single voxel: an axis ray, a face diagonal, and the body diagonal through the
   centre give `s`, `s√2`, and `s√3`.
3. Heterogeneous, to catch wrong indexing:
   - An alternating two-material slab gives `s (n_a mu_a + n_b mu_b) / 10`.
   - A 3x3x1 grid with distinct `mu_ij` and the 45° ray `b = a + 1/4`, where
     `(a, b)` are cell-unit coordinates from the box corner, gives `0.75 s√2` in
     cells (0,0), (1,1), and (2,2), `0.25 s√2` in cells (0,1) and (1,2), and
     zero elsewhere.
4. `radiological_path` equals the sum of `mu[index] * length` over `cell_path`
   bitwise, since the summation order is the same.
5. Edge cases: an exactly zero direction component; interior-face and
   upper-outer-face rays (tie rule); a vertex crossing; an inside origin
   (half-line); a miss (`None`) and a contact (`Some(0)`); non-finite `mu`,
   which propagates as documented.
6. Pose invariance: on an oriented grid with the rigidly transformed ray, the
   result equals the axis-aligned result within the bound using `eta_rot`.
7. Trapezoid identity against `forward_project_ray`. For an axis ray along a
   node line that starts outside the grid, the trilinear integral is the
   trapezoid sum, so `tau_cell - tau_march = s (v_0 + v_{n-1}) / 20`. The test
   uses a spacing `s` and a step `s / 2^k` that are exactly representable and
   asserts the plan's step count is `(n - 1) 2^k`. The midpoint plan is then
   node-aligned and exact on each affine piece.
8. Nearest-cell march: a test-local nearest-cell midpoint march over `[t0, t1]`
   satisfies `|tau_h - tau| <= (h/2) V / 10` plus roundoff, where `V` is the
   total variation of `mu` along the ray. This is asserted at several `h`, but
   not as a monotone decrease, because midpoint error on a step function
   oscillates.

## Consequences

The change is additive ([minor]), and two documented field models coexist.
Wiring the path into Radon, SIRT, or acquisition needs a consumer decision. A
GPU cell-traversal kernel belongs to Hephaestus.
