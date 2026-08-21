"""Typed interface for the ``helios`` PyO3 extension module."""

from collections.abc import Sequence


def thomson_cross_section() -> float: ...


def klein_nishina_cross_section(energy_mev: float) -> float: ...


def compton_mass_attenuation(energy_mev: float, z_over_a: float) -> float: ...


def mass_density_from_hu(hu: float, water_density_g_cm3: float) -> float: ...


def optimize_beam_weights(
    influence: Sequence[float],
    voxels: int,
    beamlets: int,
    prescription: Sequence[float],
    iterations: int,
    step: float,
) -> list[float]: ...
