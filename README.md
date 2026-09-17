\# Numerical PDEs - Course Repository

Jordan Tierney - Computational \& Applied Mathematics

Colorado School of Mines



\---





\## Overview

This repository contains all homework assignments for Numerical PDEs. 

Each assignment is fully self-contained and includes:

* MATLAB implementations
* reusable utilities
* generated figures
* full LaTeX report
* scratchwork derivations

The repository is structured so new assignments can be added easily as the course progresses.



\---



\## Repository Structure

```

numerical-pdes/
│
├── hw1/
│   ├── code/
│   ├── utils/
│   ├── figures/
│   └── report/
|
└── README.md

```

Each homework folder contains its own code, figures, utilities, and LaTeX documentation.

A small README inside each homework folder describes that specific assignment.



\---



\## Homework Index

\### \*\*HW0 - Finite Differences, ODE, Poisson (Dirichlet \& Neumann)\*\*

Folder: `HW0/`

Includes:

* first derivative
* second derivative (periodic)
* linear ODE
* Poisson equation (Dirichlet)
* Poisson equation (Neumann)
* full LaTeX report
* scratchwork derivations
* MATLAB figures and utilities

Each additional homework will follow the same structure as HW0.



\---



\## Utilities

Each homework folder contains a `utils/` directory with reusable numerical routines:

* finite-difference stencils
* ODE matrix assembly
* Poisson solvers (Dirichlet \& Neumann)
* plotting helpers

Utilities are scoped per-assignment to keep each homework self-contained.



\---



\## Running Homework Code

Navigate to the homework folder and run scripts inside `code/`:

```matlab
HW0/code/problem1\\\_first\\\_derivative
HW0/code/problem2\\\_second\\\_derivative
HW0/code/problem3\\\_ode
HW0/code/problem4\\\_poisson\\\_dirichlet
HW0/code/problem4\\\_poisson\\\_neumann
```

Figures will be saved automatically to the corresponding `figures/` folder.



\---



\## Reports and LaTeX Sources

Each homework includes:

```

report/<assignment>.pdf

report/<assignment>.tex

report/<assignment\_scratchwork>.tex

```

If the LaTeX report uses external images, they are stored in:

```

report/figures/

```

This ensures the `.tex` file can be compiled outside Overleaf.



\---



\## Author

Jordan Tierney

Computational \& Applied Mathematics

Colorado School of Mines



