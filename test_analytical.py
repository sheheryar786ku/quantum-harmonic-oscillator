import numpy as np
import pytest

from qho.analytical import exact_energies, exact_wavefunction


def test_exact_energy_levels():
    np.testing.assert_allclose(exact_energies(4), [0.5, 1.5, 2.5, 3.5])


def test_ground_state_wavefunction_formula():
    x = np.array([0.0, 1.0])
    expected = np.pi ** (-0.25) * np.exp(-x**2 / 2)
    np.testing.assert_allclose(exact_wavefunction(0, x), expected)


def test_exact_wavefunction_is_normalized_on_wide_grid():
    x = np.linspace(-8, 8, 4001)
    dx = x[1] - x[0]
    psi = exact_wavefunction(2, x)
    assert np.sum(psi**2) * dx == pytest.approx(1.0, abs=1e-6)


def test_negative_quantum_number_is_rejected():
    with pytest.raises(ValueError):
        exact_wavefunction(-1, [0.0])
