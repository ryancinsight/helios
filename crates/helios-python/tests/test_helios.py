"""Value-semantic equivalence tests for the Helios PyO3 surface.

Each assertion checks a computed value against an independent analytical oracle
(not mere importability), mirroring the Rust-side tests so the Python-visible
result is proven equal to the Rust core. Build first with `maturin develop`
in this crate (``crates/helios-python``), then ``pytest``.
"""

import math

import pytest

import helios


def test_thomson_cross_section() -> None:
    # sigma_T = 8*pi/3 * r_e^2, r_e = 2.8179403262e-15 m -> 6.652e-29 m^2.
    r_e = 2.8179403262e-15
    expected = 8.0 * math.pi / 3.0 * r_e * r_e
    assert helios.thomson_cross_section() == pytest.approx(expected, rel=1e-9)


def test_klein_nishina_thomson_limit() -> None:
    # As E -> 0, sigma_KN -> sigma_T (Thomson limit): agree to <1% at 1 keV.
    sigma_t = helios.thomson_cross_section()
    sigma_kn = helios.klein_nishina_cross_section(1e-3)
    assert sigma_kn == pytest.approx(sigma_t, rel=1e-2)


def test_klein_nishina_monotone_decreasing() -> None:
    # In the MV regime the total cross-section falls with photon energy.
    assert helios.klein_nishina_cross_section(10.0) < helios.klein_nishina_cross_section(1.0)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_klein_nishina_rejects_nonpositive(bad: float) -> None:
    with pytest.raises(ValueError):
        helios.klein_nishina_cross_section(bad)


def test_compton_mass_attenuation_water_1mev() -> None:
    # Water <Z/A> ~ 0.5551; at 1 MeV the first-principles Compton mu/rho matches
    # the NIST total (0.0707 cm^2/g), which is Compton-dominated (~99.8%) there.
    mu_rho = helios.compton_mass_attenuation(1.0, 0.5551)
    assert mu_rho == pytest.approx(0.0707, rel=2e-2)


def test_compton_mass_attenuation_rejects_nonpositive_energy() -> None:
    with pytest.raises(ValueError):
        helios.compton_mass_attenuation(-1.0, 0.5551)


def test_mass_density_from_hu_calibration() -> None:
    # Bilinear CT calibration: HU=0 -> water density, HU=1000 -> 2x water.
    assert helios.mass_density_from_hu(0.0, 1.0) == pytest.approx(1.0, abs=1e-12)
    assert helios.mass_density_from_hu(1000.0, 1.0) == pytest.approx(2.0, abs=1e-12)


def test_optimize_beam_weights_identity() -> None:
    # A = I (2x2): the unconstrained optimum equals the prescription.
    influence = [1.0, 0.0, 0.0, 1.0]
    prescription = [2.0, 3.0]
    w = helios.optimize_beam_weights(influence, 2, 2, prescription, 500, 0.5)
    assert list(w) == pytest.approx([2.0, 3.0], rel=1e-3)


def test_optimize_beam_weights_nonnegativity() -> None:
    # A negative prescription is clamped away by the x >= 0 projection.
    w = helios.optimize_beam_weights([1.0], 1, 1, [-5.0], 100, 0.5)
    assert w[0] == pytest.approx(0.0, abs=1e-9)


def test_optimize_beam_weights_rejects_shape_mismatch() -> None:
    # influence length (2) != voxels*beamlets (4).
    with pytest.raises(ValueError):
        helios.optimize_beam_weights([1.0, 2.0], 2, 2, [1.0, 1.0], 10, 0.1)
    # prescription length (1) != voxels (2).
    with pytest.raises(ValueError):
        helios.optimize_beam_weights([1.0, 0.0, 0.0, 1.0], 2, 2, [1.0], 10, 0.1)
