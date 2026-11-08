"""Command-line interface for the quantum harmonic oscillator solver."""

import argparse

import numpy as np

from .analytical import exact_energies
from .plotting import make_plots
from .solver import solve_oscillator


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Numerically solve the 1D quantum harmonic oscillator."
    )
    parser.add_argument("--states", type=int, default=6,
                        help="number of lowest states (default: 6)")
    parser.add_argument("--points", type=int, default=1000,
                        help="number of interior grid points (default: 1000)")
    parser.add_argument("--xmax", type=float, default=8.0,
                        help="half-width of domain (default: 8)")
    parser.add_argument("--output", default="results",
                        help="directory for plot files (default: results)")
    parser.add_argument("--no-plots", action="store_true",
                        help="do not save plots")
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        result = solve_oscillator(
            states=args.states, points=args.points, xmax=args.xmax
        )
    except ValueError as error:
        parser.error(str(error))

    exact = exact_energies(args.states)
    errors = np.abs(result.energies - exact)

    print("\nDimensionless energy eigenvalues")
    print("-" * 62)
    print(f"{'n':>3} {'Numerical':>18} {'Exact':>18} {'Abs. error':>18}")
    print("-" * 62)
    for n, (numerical, analytical, error) in enumerate(
        zip(result.energies, exact, errors)
    ):
        print(f"{n:3d} {numerical:18.10f} {analytical:18.10f} {error:18.3e}")
    print("-" * 62)
    print(f"Grid points: {args.points}; dx: {result.dx:.6g}; xmax: {args.xmax}")

    if not args.no_plots:
        paths = make_plots(result, args.output)
        print("\nSaved plots:")
        for path in paths:
            print(f"  {path}")

    return 0
