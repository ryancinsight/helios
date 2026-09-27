//! Allocation contract of `optimize_beam_weights`: the four work vectors are
//! allocated once, so the heap-allocation count does not grow with the
//! iteration count.
//!
//! A counting global allocator forwards to `System` and counts `alloc` calls
//! made on the measuring thread while a window is open. `realloc` and
//! `alloc_zeroed` use the trait's default bodies, which route through `alloc`,
//! so every fresh allocation is counted. The allocator is process-wide, which is
//! why this contract lives in its own test binary.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::sync::atomic::{AtomicUsize, Ordering};

use helios_planning::{optimize_beam_weights, DoseInfluence};

struct CountingAllocator;

static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

thread_local! {
    static COUNTING: Cell<bool> = const { Cell::new(false) };
}

// SAFETY: every method forwards to `System` with the caller's arguments
// unchanged; the counter only observes calls and never touches the memory.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNTING.try_with(Cell::get).unwrap_or(false) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        }
        // SAFETY: the caller upholds `GlobalAlloc::alloc`'s contract for
        // `layout`, which is passed through unchanged.
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: `ptr` was returned by `System.alloc` (via `alloc` above) for
        // this `layout`, as the caller guarantees.
        unsafe { System.dealloc(ptr, layout) }
    }
}

#[global_allocator]
static GLOBAL: CountingAllocator = CountingAllocator;

/// Heap allocations performed by one `optimize_beam_weights` call.
fn allocations_during(
    influence: &DoseInfluence<f64>,
    prescription: &[f64],
    iterations: usize,
) -> usize {
    ALLOCATIONS.store(0, Ordering::Relaxed);
    COUNTING.with(|c| c.set(true));
    let weights = optimize_beam_weights(influence, prescription, iterations, 1.0e-3);
    COUNTING.with(|c| c.set(false));
    std::hint::black_box(weights);
    ALLOCATIONS.load(Ordering::Relaxed)
}

#[test]
fn optimizer_allocates_four_work_vectors_whatever_the_iteration_count() {
    let (voxels, beamlets) = (32, 8);
    let rows: Vec<f64> = (0..voxels * beamlets)
        .map(|i| f64::from(u32::try_from(i % 13).expect("invariant: i % 13 < 13")) * 0.05)
        .collect();
    let influence = DoseInfluence::from_rows(voxels, beamlets, rows)
        .expect("invariant: rows.len() == voxels * beamlets");
    let prescription = vec![1.0; voxels];

    // x, dose, residual, gradient: allocated once before the loop. The base
    // implementation allocated three more per iteration (4 / 31 / 3001).
    let counts: Vec<usize> = [1, 10, 1000]
        .iter()
        .map(|&n| allocations_during(&influence, &prescription, n))
        .collect();
    assert_eq!(counts, vec![4, 4, 4]);
}
