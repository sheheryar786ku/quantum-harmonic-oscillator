"""Finite-difference solver for the dimensionless 1D quantum oscillator."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh_tridiagonal


@dataclass
class SolverResult:
    """Numerical grid, energies and normalized eigenstates.

    Attributes:
        x: Interior dimensionless position grid.
        dx: Uniform grid spacing.
        energies: Lowest dimensionless energy eigenvalues.
        wavefunctions: Array shaped (number_of_points, number_of_states).
                      Each column is normalized so sum(|psi|**2) * dx = 1.
    """

    x: np.ndarray
    dx: float
    energies: np.ndarray
    wavefunctions: np.ndarray


def _validate_grid(xmax: float, points: int) -> None:
    """Validate the spatial domain and number of interior grid points."""
    if not np.isfinite(xmax) or xmax <= 0:
        raise ValueError("xmax must be a finite positive number.")
    if isinstance(points, bool) or int(points) != points or points < 3:
        raise ValueError("points must be an integer of at least 3.")


def build_hamiltonian(xmax: float = 8.0, points: int = 1000):
    """Build the symmetric tridiagonal Hamiltonian in dimensionless units.

    Returns:
        x: Interior grid points.
        diagonal: Main diagonal of the Hamiltonian.
        off_diagonal: First upper/lower diagonal.
        dx: Grid spacing.

    The endpoints at -xmax and +xmax are fixed to zero (Dirichlet conditions).
    """
    _validate_grid(xmax, points)

    full_grid = np.linspace(-xmax, xmax, int(points) + 2)
    dx = float(full_grid[1] - full_grid[0])
    x = full_grid[1:-1]

    diagonal = 1.0 / dx**2 + 0.5 * x**2
    off_diagonal = np.full(int(points) - 1, -0.5 / dx**2)

    return x, diagonal, off_diagonal, dx


def solve_oscillator(
    states: int = 6, points: int = 1000, xmax: float = 8.0
) -> SolverResult:
    """Calculate the lowest energy eigenvalues and eigenstates.

    Args:
        states: Number of lowest eigenstates to calculate.
        points: Number of interior grid points.
        xmax: Positive half-width of the computational domain.

    Returns:
        SolverResult containing the grid, energies and normalized states.
    """
    if isinstance(states, bool) or int(states) != states or states < 1:
        raise ValueError("states must be a positive integer.")
    if states > points:
        raise ValueError("states cannot exceed the number of interior points.")

    x, diagonal, off_diagonal, dx = build_hamiltonian(xmax, points)

    energies, eigenvectors = eigh_tridiagonal(
        diagonal,
        off_diagonal,
        select="i",
        select_range=(0, int(states) - 1),
        check_finite=True,
    )

    # eigh_tridiagonal returns eigenvectors as columns. Normalize with
    # the discrete integral sum(|psi|^2) * dx = 1.
    norms = np.sqrt(np.sum(np.abs(eigenvectors) ** 2, axis=0) * dx)
    wavefunctions = eigenvectors / norms[np.newaxis, :]

    # Fix arbitrary eigenvector signs for consistent plots and comparisons:
    # choose each state's largest-magnitude sample to be positive.
    for column in range(wavefunctions.shape[1]):
        pivot = np.argmax(np.abs(wavefunctions[:, column]))
        if wavefunctions[pivot, column] < 0:
            wavefunctions[:, column] *= -1

    return SolverResult(
        x=x,
        dx=dx,
        energies=energies,
        wavefunctions=wavefunctions,
    )
