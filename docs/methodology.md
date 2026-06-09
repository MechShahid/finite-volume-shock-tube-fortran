# Methodology

## Governing equations

The solver computes the one-dimensional compressible Euler equations in conservative form:

```text
∂U/∂t + ∂F(U)/∂x = 0
```

where

```text
U = [rho, rho u, E]^T
F = [rho u, rho u^2 + p, u(E + p)]^T
```

The total energy is

```text
E = p/(gamma - 1) + 0.5 rho u^2
```

## Numerical method

The solution is advanced using a finite-volume update:

```text
U_i^{n+1} = U_i^n - (dt/dx) (F_{i+1/2} - F_{i-1/2})
```

The first-order `godunov` option uses the Rusanov approximate Riemann flux. It is Godunov-type in the finite-volume sense, but it is not an exact Godunov Riemann solver:

```text
F_{i+1/2} = 0.5 (F_L + F_R) - 0.5 a_max (U_R - U_L)
```

where `a_max` is the maximum local wave speed `|u| + c`.

## MUSCL reconstruction

For the MUSCL option, left and right interface states are reconstructed using a limited slope. The Van Albada limiter is used:

```text
phi(r) = (r^2 + r)/(r^2 + 1)
```

This improves spatial resolution in smooth regions while reducing spurious oscillations near discontinuities. The time update is still forward Euler, so the compact educational solver should not be described as a fully second-order production code.

## Hybrid method

The hybrid option uses MUSCL reconstruction in smooth regions and reverts to first-order reconstruction near large pressure jumps. This is a simple shock sensor approach.

## Boundary conditions

Zero-gradient/transmissive boundary conditions are applied at both ends of the domain.

## Output

Each simulation writes a CSV file with:

```text
x, rho, u, p, E
```
