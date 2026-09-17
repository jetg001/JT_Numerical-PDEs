\# Numerical PDEs - Homework 0

Jordan Tierney - Computational \& Applied Mathematics

Colorado School of Mines



\---



\## Overview

Homework 0 covers:

* First derivative approximation
* Second derivative with periodic boundary conditions
* Linear ODE with variable coefficients
* Poisson equation (Dirichlet)
* Poisson Equation (Neumann)

This folder contains all the code, figures, utilities, and LaTeX documentation for the assignment.



\---



\## Folder Structure

```

HW0/

│

├── code/

│   ├── problem1\_first\_derivative.m

│   ├── problem2\_second\_derivative.m

│   ├── problem3\_ode.m

│   ├── problem4\_poisson\_dirichlet.m

│   └── problem4\_poisson\_neumann.m

│

├── utils/

│   ├── first\_derivative\_fd.m

│   ├── second\_derivative\_periodic.m

│   ├── ode\_matrix.m

│   ├── poisson\_dirichlet\_system.m

│   ├── poisson\_neumann\_system.m

│   └── format\_loglog\_axes.m

│

├── figures/

│   ├── p1\_derivative.png

│   ├── p1\_convergence.png

│   ├── p2\_second\_derivative.png

│   ├── p2\_convergence.png

│   ├── p3\_solution.png

│   ├── p3\_convergence.png

│   ├── p4\_dirichlet\_solution.png

│   ├── p4\_dirichlet\_convergence.png

│   ├── p4\_neumann\_solution.png

│   └── p4\_neumann\_convergence.png

│

└── report/

&#x20;   ├── Numerical\_PDEs\_HW0.pdf

&#x20;   ├── Numerical\_PDEs\_HW0.tex

&#x20;   └── Numerical\_PDEs\_HW0\_scratchwork.tex

&#x20;   └── figures/

&#x09;├── p1\_derivative.png

&#x09;├── p1\_convergence.png

&#x09;├── p2\_second\_derivative.png

&#x09;├── p2\_convergence.png

&#x09;├── p3\_solution.png

&#x09;├── p3\_convergence.png

&#x09;├── p4\_dirichlet\_solution.png

&#x09;├── p4\_dirichlet\_convergence.png

&#x09;├── p4\_neumann\_solution.png

&#x09;└── p4\_neumann\_convergence.png

```



\---



\## Running the Code

Run any script inside `code/`:

```matlab

HW0/code/problem1\_first\_derivative

HW0/code/problem2\_second\_derivative

HW0/code/problem3\_ode

HW0/code/problem4\_poisson\_dirichlet

HW0/code/problem4\_poisson\_neumann

```

Figures will be saved automatically to `figures/`.



\---



\## Report \& Scratch Work

The full write-up and derivations are located in:

```

report/Numerical\_PDEs\_HW0.pdf

report/Numerical\_PDEs\_HW0.tex

report/Numerical\_PDEs\_HW0\_scratchwork.tex

```



\---



\## Notes

Homework 0 establishes a finite-difference framework for use in later assignments.

