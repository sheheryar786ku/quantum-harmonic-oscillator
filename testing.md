# Testing and validation

Run all tests from the repository root:

```bash
python -m pytest
```

The test suite verifies:

- Hamiltonian dimensions and finite-difference coefficients.
- Numerical energy eigenvalues against the exact \(n+1/2\) spectrum.
- Discrete normalization of numerical eigenstates.
- Orthogonality of the eigenstates.
- Improved ground-state accuracy after grid refinement.
- Validation of invalid parameters.
- Exact analytical energy and ground-state wavefunction formulas.
- Normalization of an analytical wavefunction on a sufficiently wide grid.

Tests are not a proof that every possible input is correct. They are repeatable
checks of important physical and numerical properties.
