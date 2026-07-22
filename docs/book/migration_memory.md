# Chapter 31 — Mnemosyne and Themis: Memory

Helios migrates from `std::alloc` (jemalloc/mimalloc on production) to
**Mnemosyne** arenas plus **Themis** NUMA placement.  Two crates compose:
Mnemosyne reuses memory in type-erased arenas; Themis binds the arenas to
the physical / logical NUMA node so that a solver thread reads its state
cache-hot.

## Mnemosyne Arenas

```rust
pub struct Arena {
    chunks: Mutex<Vec<ArenaChunk>>,
    free_list: Mutex<Vec<NonNull<u8>>>,
}

impl Arena {
    pub fn with_capacity(bytes: usize) -> Self;
    pub fn alloc<T>(&self, value: T) -> &mut T;
    pub fn try_alloc<T>(&self, value: T) -> Option<&mut T>;
    pub fn reset(&mut self);
}
```

The default helios dose-engine workflow:

1. **Allocate once at construction.** `arena: Arena<AttenuationMap>`
   instead of `Vec<f64>` outliving the engine.
2. **Reset between patients.** `arena.reset()` after each treatment plan,
   freeing everything in O(1) without per-element `drop`.
3. **Sub-arenas for transient projections.** A `ScratchArena` holds the
   per-angle intermediate `Projection`; dropped at scope exit — no `drop`
   per projection.

The helios dose engine that previously allocated millions of small
`AttenuationCoeff` structs per beam angle sees **zero heap fragmentation**
and a **predictable cache footprint** after the port.

## Themis NUMA Placement

```rust
pub struct PhysicalCore(pub u32);
pub struct NumaNode(pub u32);

pub trait Placement {
    fn bind_to(core: PhysicalCore) -> Self;
    fn numa_aware() -> Self;
    fn current() -> Self;
}
```

`themis::Placement` exposes the host topology (cores, NUMA nodes, LLC
partitions).  Helios binds its dose pools to NUMA-aware locations:

```rust
let placement = Placement::numa_aware();
let pool = MoiraiPool::new(placement.clone(), num_workers);
let arena = Arena::with_capacity(dose_grid_bytes).bind(placement);
```

The binding does not move memory at runtime — it influences which **arena
chunk** the allocator serves next, so consecutive allocations land on
the same NUMA node and the solver reads cache-hot data.

## Migration Procedure

| Legacy | Atlas |
|---|---|
| global `jemalloc`/`mimalloc` init | `Arena::with_capacity(...)` per subsystem |
| `Vec<T>::with_capacity(N)` | `Arena::alloc_slice::<T>(N)` |
| per-step `drop` cost | `arena.reset()` once per patient |
| `sched_setaffinity` (Linux) | `Placement::bind_to(core)` (portable) |

## Validation Examples

- [`tomotherapy_workflow`](examples/tomotherapy_workflow.md) —
  per-patient arena lifecycle for thousands of beams.
- [`adaptive_rt_workflow`](examples/adaptive_rt_workflow.md) — arena
  reuse across adaptive replanning cycles.
- [`linac_dose_accumulation`](examples/linac_dose_accumulation.md) —
  MLC-aware arena pool.
- [`gpu_attenuation_projection`](examples/gpu_attenuation_projection.md) —
  host-side scratch arena in the GPU forward path.

## Further Reading

- [`mnemosyne` source](../../../mnemosyne/crates/)
- [`themis` source](../../../themis/src/)
- [`consus` source](../../../consus/) — persistent storage on top of these.
