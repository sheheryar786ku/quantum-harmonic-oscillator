import numpy as np
import pytest

from qho.solver import build_hamiltonian, solve_oscillator


def test_hamiltonian_diagonal_and_off_diagonal():
    points = 10
    xmax = 4.0
    x, diagonal, off_diagonal, dx = build_hamiltonian(xmax, points)

    assert len(x) == points
    assert len(diagonal) == points
    assert len(off_diagonal) == points - 1
    np.testing.assert_allclose(diagonal, 1.0 / dx**2 + 0.5 * x**2)
    np.testing.assert_allclose(off_diagonal, -0.5 / dx**2)


def test_lowest_energies_match_exact_spectrum():
    result = solve_oscillator(states=4, points=1200, xmax=8.0)
    exact = np.arange(4) + 0.5
    # The finite-difference error at this resolution is comfortably below 0.01.
    np.testing.assert_allclose(result.energies, exact, atol=0.002, rtol=0)


def test_wavefunctions_are_discretely_normalized():
    result = solve_oscillator(states=5, points=500, xmax=8.0)
    norms = np.sum(result.wavefunctions**2, axis=0) * result.dx
    np.testing.assert_allclose(norms, np.ones(5), atol=1e-12)


def test_eigenstates_are_orthogonal():
    result = solve_oscillator(states=5, points=500, xmax=8.0)
    overlap = result.wavefunctions.T @ result.wavefunctions * result.dx
    np.testing.assert_allclose(overlap, np.eye(5), atol=1e-11)


def test_grid_refinement_improves_ground_state_energy():
    coarse = solve_oscillator(states=1, points=100, xmax=8.0).energies[0]
    fine = solve_oscillator(states=1, points=400, xmax=8.0).energies[0]
    exact = 0.5
    assert abs(fine - exact) < abs(coarse - exact)


@pytest.mark.parametrize("xmax,points", [(0, 10), (-1, 10), (3, 2)])
def test_invalid_grid_is_rejected(xmax, points):
    with pytest.raises(ValueError):
        build_hamiltonian(xmax=xmax, points=points)


def test_invalid_state_count_is_rejected():
    with pytest.raises(ValueError):
        solve_oscillator(states=0)
