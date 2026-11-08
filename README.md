# Numerical Quantum Harmonic Oscillator Solver

A documented, tested Python project for the **one-dimensional quantum harmonic
oscillator**. It approximates the time-independent Schrödinger equation using a
finite-difference Hamiltonian, calculates energy eigenvalues and eigenstates,
and compares numerical energies with the exact analytical spectrum.

This project is designed as a computational-physics programming project. It
includes a command-line interface, reusable source modules, automated tests,
plots, mathematical background, and reproducible installation instructions.

## Physics model

The time-independent Schrödinger equation is

\[
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}
+\frac12m\omega^2x^2\psi=E\psi.
\]

We use dimensionless variables

\[
q=x\sqrt{\frac{m\omega}{\hbar}},\qquad
\epsilon=\frac{E}{\hbar\omega},
\]

so that the equation becomes

\[
\left[-\frac12\frac{d^2}{dq^2}+\frac12q^2\right]\psi(q)
=\epsilon\psi(q).
\]

The exact dimensionless energies are \(\epsilon_n=n+1/2\), where
\(n=0,1,2,\ldots\).

On a uniform grid with spacing \(h\), the second derivative is approximated by

\[
\psi''(q_i)\approx
\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{h^2}.
\]

This produces a real symmetric tridiagonal Hamiltonian. Its diagonal entries
are \(1/h^2+q_i^2/2\), and its two adjacent off-diagonals are
\(-1/(2h^2)\). We diagonalize this matrix with
`scipy.linalg.eigh_tridiagonal`.

The numerical problem is solved on a finite interval with Dirichlet boundary
conditions (the wavefunction is zero at the two endpoints). Results therefore
depend slightly on grid spacing and interval size. Increase `--points` and
`--xmax` to check convergence.

## Features

- Finite-difference Hamiltonian construction.
- Efficient symmetric tridiagonal eigensolver.
- Numerical energy eigenvalues and eigenstates.
- Discrete wavefunction normalization.
- Exact analytical energies and Hermite-function wavefunctions for comparison.
- Energy-error table and plots of energy levels, wavefunctions and probability densities.
- Command-line options for grid size, domain and number of states.
- Automated unit tests for matrix structure, energy accuracy, normalization,
  orthogonality and convergence.

## Requirements

- Python 3.9 or newer
- NumPy
- SciPy
- Matplotlib
- pytest (for running tests)

## Installation

Clone the public repository (replace the URL with your own GitHub URL):

```bash
git clone https://github.com/YOUR-USERNAME/quantum-harmonic-oscillator.git
cd quantum-harmonic-oscillator
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate       # macOS/Linux
# .venv\Scripts\activate        # Windows
```

Install the project and development dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Alternatively, install from `requirements.txt`:

```bash
python -m pip install -r requirements.txt
```

## Run the solver

```bash
python -m qho --states 6 --points 1000 --xmax 8
```

The command prints a table comparing numerical energies with exact values
and saves three plots in the `results/` directory.

You can also use the installed command:

```bash
qho-solver --states 6 --points 1000 --xmax 8
```

### Command-line options

| Option | Default | Meaning |
|---|---:|---|
| `--states` | 6 | Number of lowest states to calculate |
| `--points` | 1000 | Number of interior grid points |
| `--xmax` | 8.0 | Domain is approximately `[-xmax, xmax]` |
| `--output` | results | Directory for generated plots |
| `--no-plots` | off | Print energies without making plots |

Example for a higher-resolution calculation:

```bash
python -m qho --states 8 --points 1600 --xmax 9
```

## Run the tests

From the repository root:

```bash
python -m pytest
```

The tests check:
1. Hamiltonian diagonal and off-diagonal elements.
2. Agreement of the first few numerical energies with the exact spectrum.
3. Discrete normalization of each computed eigenstate.
4. Orthogonality of different eigenstates.
5. Improvement of the ground-state energy as the grid is refined.

## Project structure

```text
quantum-harmonic-oscillator/
├── README.md
├── LICENSE
├── pyproject.toml
├── requirements.txt
├── .gitignore
├── src/
│   └── qho/
│       ├── __init__.py
│       ├── __main__.py
│       ├── solver.py
│       ├── analytical.py
│       ├── plotting.py
│       └── cli.py
├── tests/
│   ├── test_solver.py
│   └── test_analytical.py
├── docs/
│   ├── theory.md
│   ├── numerical_method.md
│   └── testing.md
└── results/
    └── .gitkeep
```

## Reproducibility and limitations

The calculation uses a finite spatial interval and a finite grid. This introduces
domain truncation and discretization errors. The reported energies are
dimensionless, in units of \(\hbar\omega\). To convert to physical energy,
multiply by \(\hbar\omega\). The dimensionless coordinate is
\(q=x\sqrt{m\omega/\hbar}\).

The discrete normalization convention is
\(\sum_i |\psi_i|^2 h=1\). The returned wavefunctions are normalized using
this convention and have arbitrary overall signs, as expected for eigenvectors.

## Suggested future extensions

- Interactive state selection and parameter controls.
- Convergence studies over grid spacing and domain size.
- Anharmonic potentials such as \(V(q)=q^2/2+\lambda q^4\).
- Time-dependent wavepacket evolution.
- Continuous-integration testing with GitHub Actions.

## Academic use

Use this repository as a starting point, but make sure you can explain the
physics, finite-difference approximation, matrix construction, numerical
errors, tests and code design during your examination. Follow your course's
rules about attribution and permitted collaboration.
