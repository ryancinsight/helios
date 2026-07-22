# Appendix A — Atlas Crate Dependency Map

Helios is an Atlas domain consumer. Dependency direction runs from workflow
orchestration toward domain vocabulary and provider-owned infrastructure.

Primary production edges in the workspace manifests are:

- `helios-simulation` → `helios-solver`, `helios-physics`, and
  `helios-domain`;
- `helios-imaging` → `helios-solver` and `helios-domain`;
- `helios-solver` → `helios-physics` and `helios-domain`;
- `helios-analysis` → `helios-domain`; and
- `helios-domain` → `helios-math` and `helios-core`.

The end-to-end examples and tests additionally compose simulation, imaging, and
analysis through development-only dependencies. Atlas providers own numeric
traits (`eunomia`), arrays and linear algebra (`leto`), spatial primitives
(`gaia`), photon transport (`hyperion`), GPU compute (`hephaestus`), parallel
execution (`moirai`), DICOM decoding (`ritk-dicom`), and volume persistence
(`consus-hdf5`).

The root [architecture document](../../ARCHITECTURE.md) owns the complete
layering and dependency contract. This appendix is a navigation aid, not a
second dependency specification.

## Further reading

- [API Reference Index](appendix_api.md)
