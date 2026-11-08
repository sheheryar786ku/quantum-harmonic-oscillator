# Numerical method

## 1. Spatial grid

The finite interval \([-q_{\max},q_{\max}]\) is divided into `points + 1`
equal intervals. The endpoints are excluded from the unknown vector and the
wavefunction is set to zero at both boundaries (Dirichlet boundary conditions).

## 2. Finite-difference approximation

At an interior grid point \(q_i\),

\[
\psi''(q_i)\approx
\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{\Delta q^2}.
\]

Substituting this into the dimensionless Hamiltonian gives

\[
H_{ii}=\frac{1}{\Delta q^2}+\frac{q_i^2}{2},
\qquad
H_{i,i+1}=H_{i+1,i}=-\frac{1}{2\Delta q^2}.
\]

All other matrix elements are zero. The Hamiltonian is therefore real,
symmetric and tridiagonal.

## 3. Diagonalization

The eigenvalue problem is \(H\mathbf{v}_n=\epsilon_n\mathbf{v}_n\).
`scipy.linalg.eigh_tridiagonal` computes the lowest requested eigenpairs
without constructing a dense matrix. This is more memory-efficient for a
tridiagonal problem.

The eigenvectors are rescaled so that
\(\sum_i|v_{i,n}|^2\Delta q=1\), matching the continuum normalization integral.

## 4. Sources of numerical error

- **Grid discretization:** the finite-difference formula has leading error
  proportional to \((\Delta q)^2\).
- **Finite domain:** the exact wavefunction extends to infinity, while the
  computation imposes zero values at finite endpoints.
- **Floating-point arithmetic:** calculations use finite-precision numbers.

A convergence study should vary grid resolution and domain size independently.
For low-lying states, a sufficiently large domain and finer grid should bring
the numerical energies closer to \(n+1/2\).
