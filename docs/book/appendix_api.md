# Appendix C — API Reference Index

The Rustdoc output is the authoritative API reference:

```bash
cargo doc --workspace --no-deps --open
```

## Selected public APIs

### `helios-core`

- `EnergyMeV`, `HounsfieldUnit`, and `VoxelSpacingMm`
- `HeliosError` and `Result`

### `helios-domain`

- `VoxelGrid<T>` and `Volume<T>`
- `HelicalDelivery<T>`, `LeafOpenTimeSinogram<T>`, and `MlcModel<T>`
- `FieldAperture<T>`

### `helios-imaging`

- `parallel_beam_radon` and `filtered_back_projection`
- `sirt_reconstruction`
- `register_translation` and `register_translation_ncc`

### `helios-simulation`

- `simulate_helical_delivery` and `simulate_helical_sinogram`
- `accumulate_delivered_dose` and `accumulate_delivered_dose_anisotropic`
- `BeamGeometry` and `CollapsedCone`

### `helios-analysis`

- `Dvh`
- `gamma_index_3d`, `gamma_index_3d_local`, and `gamma_pass_rate`
- Raw MVCT metrics: `roi_statistics`, `volume_rmse`, and
  `volume_relative_l2_error`
- Dose-semantic metrics: `dose_roi_statistics` and `dose_volume_rmse`

## Further reading

- [Repository changelog](../../CHANGELOG.md)
