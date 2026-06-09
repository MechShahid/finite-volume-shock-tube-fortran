\# Finite-Volume Simulation of the 1D Sod Shock Tube Problem in Fortran



This repository contains a cleaned Fortran implementation of the classical one-dimensional Sod shock tube problem for the compressible Euler equations.



The original assignment compared velocity, density, and pressure profiles for mesh sizes `nx = 201` and `nx = 1001` at `t = 0.25 s` using Godunov, MUSCL, hybrid Godunov-MUSCL, and exact solutions. This repository keeps the original assignment material for reference and provides a cleaner, reproducible, GitHub-ready version of the solver.



\## Methods Included



\* `godunov`: first-order Godunov-type finite-volume method using the Rusanov approximate Riemann flux

\* `muscl`: MUSCL reconstruction with the Van Albada limiter

\* `hybrid`: MUSCL reconstruction in smooth regions and first-order flux near strong pressure gradients



\## Problem Setup



The computational domain is:



```text

x in \[0, 1]

```



The Sod shock tube initial condition is:



```text

Left state:  rho = 1.0,   u = 0.0, p = 1.0

Right state: rho = 0.125, u = 0.0, p = 0.1



gamma = 1.4

final time = 0.25 s

```



\## Repository Structure



src/              Fortran solver

scripts/          Python plotting and exact-solution comparison script

results/          Generated CSV result files

figures/          Generated comparison figures

report/           Technical report and supporting documentation



\## Requirements



\* `gfortran`

\* `make`

\* Python 3

\* NumPy

\* Matplotlib



\## Build



Compile using:



```bash

make

```



Or compile manually:



```bash

gfortran -O2 -Wall -Wextra -std=f2008 src/shock\_tube\_solver.f90 -o shocktube

```



\## Run All Standard Cases



```bash

make run

```



This generates results for:



```text

nx = 201 and 1001

schemes = godunov, muscl, hybrid

final time = 0.25 s

```



\## Run One Case Manually



```bash

./shocktube godunov 201 0.25 results/godunov\_n201.csv

./shocktube muscl   201 0.25 results/muscl\_n201.csv

./shocktube hybrid  201 0.25 results/hybrid\_n201.csv

```



General command format:



```bash

./shocktube <scheme> <nx> <final\_time> <output\_file>

```



\## Plot Comparisons with the Exact Solution



```bash

make plots

```



The comparison figures are written to:



```text

figures/

```



\## Implementation Features



The cleaned GitHub version includes:



\* a single command-line Fortran solver with selectable numerical scheme

\* adaptive CFL-based time stepping

\* robust MUSCL reconstruction with Van Albada limiter

\* conservative-variable-based density handling

\* positivity safeguards for density and pressure

\* CSV output files for post-processing and plotting

\* Python script for exact-solution comparison

\* Makefile for reproducible compilation and execution

\* technical report and methodology documentation



\## Author



Mohammad Shahid

##### 

