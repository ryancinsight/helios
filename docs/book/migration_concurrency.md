# Moirai: Concurrency — Executors, Task Graphs, Channels

Helios migrates parallel execution from `rayon::ParIter`, bespoke
`std::thread::spawn`, and partial `tokio` usage to **Moirai**'s unified
executor.  Moirai unifies what tokio does for async (DICOM streaming)
and what rayon does for parallel numerical kernels (FBP, MVCT, DVH), so
helios workloads share one runtime.

## Executor Surface

```rust
pub struct Executor {
    inner: Arc<ExecutorInner>,
}

impl Executor {
    pub fn new(num_workers: usize) -> Self;
    pub fn num_cpus() -> usize;
    pub fn block_on<F: Future>(&self, fut: F) -> F::Output;
    pub fn spawn<F: Future + Send + 'static>(&self, fut: F) -> JoinHandle<F::Output>;
}
```

The runtime presents an `Executor` whose `.spawn` and `.block_on` accept
**any** future — async I/O and parallel numerical kernels share the
same worker pool.

## Task Graph

```rust
pub struct TaskGraph {
    nodes: Vec<TaskNode>,
    edges: Vec<(TaskId, TaskId)>,
}

impl TaskGraph {
    pub fn add_task(&mut self, ...) -> TaskId;
    pub fn add_edge(&mut self, from: TaskId, to: TaskId);
    pub fn execute(&self, exec: &Executor) -> Vec<TaskResult>;
}
```

`TaskGraph` is the migration surface for helios when a workflow has
explicit dependencies (e.g. forward-project → solve → DVH → DVH-constrained
optimization).

## Channels

```rust
pub fn channel<T: Send>(cap: usize) -> (Sender<T>, Receiver<T>);
```

Moirai channels are MPMC and back-pressure-aware.  Helios uses them for
streaming DICOM frames into the reconstruction engine without blocking
the executor.

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `rayon::par_iter().for_each(...)` | `exec.scope(\|s\| s.spawn_many(...))` |
| `rayon::join(a, b)` | `TaskGraph::add_edge` + `execute` |
| `std::thread::spawn` | `Executor::spawn` |
| `tokio::spawn` | `Executor::spawn` |
| `tokio::sync::mpsc::channel` | `moirai::channel<T>(cap)` |
| `std::sync::mpsc` | `moirai::channel<T>(cap)` |

A typical helios parallel region becomes:

```rust
exec.scope(|scope| {
    for beam in plan.beams() {
        scope.spawn(async move { compute_dose(beam, &state) });
    }
});
```

The executor decides worker count, schedules across cores, and avoids
oversubscription — three wins the legacy code reproduced per-beam.

## Validation Examples

- [`tomotherapy_workflow`](examples/tomotherapy_workflow.md) —
  TomoTherapy helical reconstruction parallelized via `TaskGraph`.
- [`linac_dose_accumulation`](examples/linac_dose_accumulation.md) —
  per-beam parallel dose accumulation.
- [`adaptive_rt_workflow`](examples/adaptive_rt_workflow.md) — adaptive
  replanning, parallel-over-runs.
- [`sirt_reconstruction`](examples/sirt_reconstruction.md) — iterative
  reconstruction under Moirai executors.

## Further Reading

- [`moirai` source](../../../moirai/crates/)
- [Mnemosyne and Themis: Memory](migration_memory.md)
- [Moirai Book cross-reference](../../../moirai/docs/)
