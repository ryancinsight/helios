# Chapter 34 — Leto: GAT-Based Tiling and Lending Iterators

Helios iterates over volumetric data in **tiles** (sliding windows,
rotational/projection haloes, beam-rotation carousels) at multiple
granularities.  The legacy approach allocates a `Vec<Tile>` per call and
pays the heap traffic.  Atlas **Leto** uses **generic associated types
(GATs)** to encode tile iteration without an internal allocator.

## The GAT Tile Iterator

```rust
pub trait TileStreaming<'a> {
    type Item;
    type LendingIter: LendingIterator<Item = &'a Self::Item> + 'a;
    fn tiles(&'a self) -> Self::LendingIter;
}
```

`TileStreaming::tiles` returns a **lending iterator** that borrows from
the source.  Each `next()` yields a `&Tile<'a>` that lives only as long as
the iterator state — the type system enforces no `'static` clones, so the
hot path is **zero-copy and zero-allocator**.

## The LendingIterator Trait

```rust
pub trait LendingIterator {
    type Item<'a> where Self: 'a;
    fn next(&mut self) -> Option<Self::Item<'_>>;
}
```

`LendingIterator::Item` is a **GAT** — its lifetime parameter lets the
iterator yield references whose lifetime is **shorter** than `&self`,
which is impossible with `Iterator` (whose `Item` is fixed at the trait
level).

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `Vec<Tile>` allocation per pass | `TileStreaming::tiles(&dose_grid)` |
| `impl Iterator<Item = &Tile<'a>>` | `LendingIterator<Item<'a> = &Tile<'a>>` |
| per-tile `clone()` of the inner buffer | `&'a Tile<'a>` borrow |
| ghost-cell halo allocation | `halo_lazy(window_offset, Self::slice)` |

A typical helios port:

```rust
for tile in dose_grid.tiles() {
    // tile: &Tile<'_>, borrowed from `dose_grid`
    update_dose(tile, beam)?;
}
// No tile allocated, no per-step clone.
```

## How Hermes Composes

When `TileStreaming::tiles` is iterated over an `NdArray<F, Ix3>`, the
tile inner slices are `&[F]` — flat enough for [`hermes-simd`] to apply
vectorized projections without per-call conversion.  The Leto + Hermes
combination keeps the kernel-loop short and the SIMD throughput high.

## Validation Examples

- [`dvh_analysis`](examples/dvh_analysis.md) — histogram tiles per
  region.
- [`gamma_index`](examples/gamma_index.md) — sliding-window comparison
  tiles.
- [`tomotherapy_workflow`](examples/tomotherapy_workflow.md) — per-rotation
  dose tiles.
- [`gpu_attenuation_projection`](examples/gpu_attenuation_projection.md) —
  tile-streamed projection queue.

## Further Reading

- [`leto` GAT module](../../../leto/crates/)
- [Leto: Arrays](migration_arrays.md)
- [Migration Overview](migration_overview.md)
