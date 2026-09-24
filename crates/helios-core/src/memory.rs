//! Memory allocation: the Atlas memory subsystem, adopted once for the workspace.
//!
//! The Atlas stack ships its own user-space allocator, Mnemosyne, and Helios
//! adopts it as the workspace allocation package rather than leaving every crate
//! to reach for the `std` default. This module is the seam. It is the only place
//! in Helios that names the `mnemosyne` package, so the allocator can be changed
//! by editing one file, and higher layers depend on [`helios_core::memory`](crate::memory)
//! instead of on `mnemosyne` directly.
//!
//! Two groups are re-exported:
//!
//! - `Mnemosyne` and `MnemosyneAllocator` — the `GlobalAlloc` implementations.
//!   A *binary* adopts one once with [`install_global_allocator!`]; a library
//!   never installs it, because choosing the allocator is the application's
//!   decision rather than a dependency's.
//! - The scratch surface — `ScratchPool`, `ScratchBank`, `AlignedVec` —
//!   thread-local reusable buffers for hot paths that would otherwise allocate
//!   per call.
//!
//! # Why a global allocator rather than an allocator-backed container
//!
//! `leto` (the array substrate Helios builds its volumes on) exposes a
//! `mnemosyne-alloc` feature that backs `Array` with `MnemosyneStorage`. Helios
//! does not use it, for two reasons. `MnemosyneStorage` aligns to
//! `align_of::<T>()`, exactly as `Vec` does, so it buys no cache-line property
//! that `AlignedVec` does not already provide; and it implements neither `Clone`
//! nor `Debug`, while `Volume` derives both, so adopting it would remove a
//! public capability to obtain a property a global allocator already delivers.
//! Installing the global allocator routes *every* allocation in the program —
//! `Vec`, `Box`, and `VecStorage` inside `leto::Array` alike — through Mnemosyne,
//! which is a strictly larger claim than swapping one container's storage.
//!
//! # Measuring
//!
//! Call `warm_current_thread` before opening a measurement window. It flushes
//! thread-local initialisation traffic — options parsing, arena segment
//! acquisition, per-thread allocator setup — so a timing run measures the code
//! under test rather than the allocator's first-touch cost.
//!
//! This module is gated on the default-on `mnemosyne-memory` feature.

pub use mnemosyne::{
    warm_current_thread, AlignedBuf, AlignedVec, Mnemosyne, MnemosyneAllocator, ScratchBank,
    ScratchElement, ScratchPool, DEFAULT_SCRATCH_ALIGN,
};

/// Installs [`Mnemosyne`] as the current binary's global allocator.
///
/// Expands to a single `#[global_allocator]` static, so it must be invoked at
/// most once per binary — at the crate root of a `main.rs`, an example, a
/// benchmark, or a `cdylib` — and never from a library, which would impose the
/// choice on every binary downstream of it.
///
/// ```rust,ignore
/// helios_core::install_global_allocator!();
///
/// fn main() { /* every allocation now goes through Mnemosyne */ }
/// ```
#[macro_export]
macro_rules! install_global_allocator {
    () => {
        #[global_allocator]
        static HELIOS_GLOBAL_ALLOCATOR: $crate::memory::Mnemosyne = $crate::memory::Mnemosyne;
    };
}

#[cfg(test)]
mod tests {
    use super::{ScratchPool, DEFAULT_SCRATCH_ALIGN};

    #[test]
    fn scratch_pool_hands_out_an_aligned_exact_length_buffer() {
        std::thread_local! {
            static POOL: ScratchPool<f64> = const { ScratchPool::new() };
        }

        POOL.with(|pool| {
            pool.with_scratch(37, |scratch| {
                assert_eq!(scratch.len(), 37);
                let base = scratch.as_ptr() as usize;
                assert_eq!(
                    base % DEFAULT_SCRATCH_ALIGN,
                    0,
                    "scratch base {base:#x} is not {DEFAULT_SCRATCH_ALIGN}-byte aligned",
                );
                scratch.fill(1.0);
            });
        });
    }

    #[test]
    fn warm_current_thread_is_idempotent() {
        super::warm_current_thread();
        super::warm_current_thread();
    }
}
