"""Plotting helpers for solver output."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .analytical import exact_wavefunction
from .solver import SolverResult


def make_plots(result: SolverResult, output_dir="results") -> list[Path]:
    """Save energy, wavefunction and probability-density plots."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths = []

    # Energy spectrum
    n = np.arange(len(result.energies))
    exact = n + 0.5
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(np.zeros_like(n), result.energies, label="Numerical", zorder=3)
    ax.scatter(np.ones_like(n), exact, marker="x", label="Exact", zorder=3)
    for i, energy in enumerate(result.energies):
        ax.plot([-0.08, 1.08], [energy, energy], alpha=0.2)
        ax.annotate(f"n={i}", (0.05, energy), xytext=(4, 4),
                    textcoords="offset points")
    ax.set_xticks([0, 1], ["Numerical", "Analytical"])
    ax.set_ylabel(r"Dimensionless energy $E/(\hbar\omega)$")
    ax.set_title("Quantum harmonic oscillator energy spectrum")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path = output / "energy_levels.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    # Wavefunctions
    fig, ax = plt.subplots(figsize=(9, 6))
    for state in range(len(result.energies)):
        numerical = result.wavefunctions[:, state]
        analytical = exact_wavefunction(state, result.x)
        ax.plot(result.x, numerical, label=f"Numerical n={state}")
        ax.plot(result.x, analytical, "--", linewidth=1,
                label=f"Exact n={state}")
    ax.set_xlabel("Dimensionless position")
    ax.set_ylabel(r"Wavefunction $\psi_n(x)$")
    ax.set_title("Numerical and analytical wavefunctions")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    path = output / "wavefunctions.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    # Probability densities
    fig, ax = plt.subplots(figsize=(9, 6))
    for state in range(len(result.energies)):
        probability = np.abs(result.wavefunctions[:, state]) ** 2
        ax.plot(result.x, probability, label=f"n={state}")
    ax.set_xlabel("Dimensionless position")
    ax.set_ylabel(r"Probability density $|\psi_n(x)|^2$")
    ax.set_title("Quantum probability densities")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path = output / "probability_densities.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    return paths
