//! The Helios program allocator.
//!
//! Helios programs — the examples and criterion benches — allocate through
//! Mnemosyne, the Atlas user-space allocator. A program opts in with one line at
//! its crate root:
//!
//! ```
//! helios_allocator::install_global_allocator!();
//!
//! fn main() {
//!     // Every allocation below goes through Mnemosyne.
//!     let samples: Vec<u64> = (1..=1024).collect();
//!     assert_eq!(samples.iter().sum::<u64>(), 524_800);
//! }
//! ```
//!
//! Only programs depend on this crate, and they reach it through
//! `[dev-dependencies]` (examples, benches) or `[dependencies]` (binaries).
//! Library crates never do: choosing the global allocator is the program's
//! decision, and a library dependency would make every consumer compile the
//! allocator whether it installs it or not (ADR 0018). `helios-python` is a
//! `cdylib` loaded into the Python interpreter and must not replace its
//! allocator, so it does not depend on this crate either.
//!
//! The `GlobalAlloc` contract — requested alignment honoured, null for a
//! zero-size request, `realloc` preserving contents — is tested in Mnemosyne
//! itself, at `crates/mnemosyne/tests/global_alloc_tests/{basic,policy,realloc}.rs`.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

pub use mnemosyne::Mnemosyne;

/// Installs [`Mnemosyne`] as the current program's global allocator.
///
/// Expands to one `#[global_allocator]` static, so it is invoked exactly once
/// per program, at the crate root of a binary, example, or benchmark.
#[macro_export]
macro_rules! install_global_allocator {
    () => {
        #[global_allocator]
        static HELIOS_GLOBAL_ALLOCATOR: $crate::Mnemosyne = $crate::Mnemosyne;
    };
}
