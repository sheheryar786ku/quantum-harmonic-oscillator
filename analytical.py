"""Exact dimensionless harmonic-oscillator energies and wavefunctions."""

import numpy as np
from scipy.special import eval_hermite, factorial


def exact_energies(states: int) -> np.ndarray:
    """Return exact dimensionless energies epsilon_n = n + 1/2."""
    if isinstance(states, bool) or int(states) != states or states < 1:
        raise ValueError("states must be a positive integer.")
    n = np.arange(int(states), dtype=float)
    return n + 0.5


def exact_wavefunction(n: int, x) -> np.ndarray:
    """Return normalized dimensionless eigenfunction psi_n(x).

    Uses physicists' Hermite polynomials:
        psi_n(x) = pi^(-1/4) / sqrt(2^n n!) * H_n(x) * exp(-x^2/2).
    """
    if isinstance(n, bool) or int(n) != n or n < 0:
        raise ValueError("n must be a non-negative integer.")
    x = np.asarray(x, dtype=float)
    prefactor = np.pi ** (-0.25) / np.sqrt((2.0 ** int(n)) * factorial(int(n)))
    return prefactor * eval_hermite(int(n), x) * np.exp(-0.5 * x**2)
