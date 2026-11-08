"""Numerical quantum harmonic oscillator package."""

from .solver import SolverResult, build_hamiltonian, solve_oscillator

__all__ = ["SolverResult", "build_hamiltonian", "solve_oscillator"]
__version__ = "1.0.0"
