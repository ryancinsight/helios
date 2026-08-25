# helios-python

Thin PyO3 bindings for the Helios radiation-therapy simulation and imaging
stack. The Rust crates own the physical laws, validation, and numerical
kernels; this package converts Python values at the boundary and exposes the
same value-semantic results through `import helios`.

## Install

```bash
python -m pip install helios-python
```

The wheels use the stable CPython 3.9 ABI. One wheel per supported operating
system serves CPython 3.9 and newer.

The wheel includes `helios.pyi` and the PEP 561 `py.typed` marker, so type
checkers resolve the five exported functions without importing the extension
at analysis time.

## Example

```python
from helios import klein_nishina_cross_section, optimize_beam_weights

cross_section = klein_nishina_cross_section(1.0)
weights = optimize_beam_weights(
    influence=[1.0, 0.0, 0.0, 1.0],
    voxels=2,
    beamlets=2,
    prescription=[1.0, 2.0],
    iterations=32,
    step=0.1,
)
assert cross_section > 0.0
assert len(weights) == 2
```

The Python layer is intentionally thin: compute-heavy calls release the GIL
and delegate to the typed Helios physics and planning crates. The source-backed
teaching book and runnable Rust examples are published at
<https://ryancinsight.github.io/helios/>.
