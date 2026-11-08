"""Example script: run the solver and save plots."""
from qho.solver import solve_oscillator
from qho.analytical import exact_energies
from qho.plotting import make_plots
import numpy as np

result = solve_oscillator(states=6, points=1000, xmax=8.0)
print(" n     numerical       exact       absolute error")
for n, (num, exact) in enumerate(zip(result.energies, exact_energies(6))):
    print(f"{n:2d}  {num:12.8f}  {exact:10.8f}  {abs(num-exact):.3e}")
make_plots(result, "results")
