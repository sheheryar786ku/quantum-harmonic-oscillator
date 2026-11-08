# Physical theory

## 1. Potential and Hamiltonian

For a particle of mass \(m\) in a one-dimensional harmonic potential,

\[
V(x)=\frac12m\omega^2x^2.
\]

The Hamiltonian is the sum of kinetic and potential energy:

\[
\hat H=-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}
+\frac12m\omega^2x^2.
\]

The stationary Schrödinger equation is \(\hat H\psi_n=E_n\psi_n\).

## 2. Dimensionless form

Define \(q=x\sqrt{m\omega/\hbar}\) and
\(\epsilon=E/(\hbar\omega)\). The equation becomes

\[
\left[-\frac12\frac{d^2}{dq^2}+\frac12q^2\right]\psi_n(q)
=\epsilon_n\psi_n(q).
\]

This removes the need to choose a particular particle mass or oscillator
frequency. Energies produced by the program are in units of \(\hbar\omega\).

## 3. Exact solution

The exact energy spectrum is

\[
\epsilon_n=n+\frac12,\qquad n=0,1,2,\ldots
\]

The nonzero ground-state energy \(\epsilon_0=1/2\) is called zero-point energy.

The normalized eigenfunctions are

\[
\psi_n(q)=
\frac{\pi^{-1/4}}{\sqrt{2^n n!}}
H_n(q)e^{-q^2/2},
\]

where \(H_n\) is the physicists' Hermite polynomial.

## 4. Probability interpretation

The probability density is \(|\psi_n(q)|^2\). For a normalized state,

\[
\int_{-\infty}^{\infty}|\psi_n(q)|^2\,dq=1.
\]

The code uses a finite grid, so it checks the discrete approximation
\(\sum_i|\psi_n(q_i)|^2\Delta q=1\).
