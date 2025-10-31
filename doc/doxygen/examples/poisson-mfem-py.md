@page poisson-mfem-py Mixed FE discretization of the Poisson equation (Python)
<!--
SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file AUTHORS.md
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception or LGPL-3.0-or-later
-->

The following example explains how to solve a mixed formulation of the Poisson problem
with Lagrange elements for the solution and Raviart-Thomas elements for the grid.
The full example file @ref source-poisson-mfem-py can be found in the `examples/` subdirectory.

Be warned that writing a quadrature loop for a stiffness matrix with [dune-functions][]
and Python is currently quite slow. This approach is therefore mainly useful
for small problems, for educational purposes, or if you need complete control.
For faster and more efficient ways to assemble finite element stiffness matrices,
have a look at the [dune-fufem](https://gitlab.dune-project.org/fufem/dune-fufem)
and [dune-fem](https://gitlab.dune-project.org/dune-fem/dune-fem) modules.

## Mixed formulation of the Poisson problem

Let \f$ \Omega \f$ be a bounded open subset of \f$ \mathbb{R}^2 \f$ or \f$ \mathbb{R}^3 \f$.
The Poisson problem asks for a scalar function \f$ u : \Omega \to \mathbb{R} \f$
that solves
\f[
 -\Delta u = f
 \qquad
 \text{on $\Omega$}
\f]
for a given function \f$f : \Omega \to \mathbb{R} \f$. To make the solution unique,
it is additionally required that the solution \f$ u \f$ satisfies boundary conditions
\f[
 u = g
 \qquad
 \text{on $\partial \Omega$}.
\f]

TODO: Finish this...

<!-- Include the entire file, that's better than not showing anything at all -->

@include poisson-mfem.py


[dune]: https://dune-project.org
[dune-functions]: https://gitlab.dune-project.org/staging/dune-functions


<div class="section_buttons">
| Previous      |                         Next |
|:--------------|-----------------------------:|
| @ref poisson-pq2-py |                       |
</div>



@page source-poisson-mfem-py poisson-mfem.py

This is the raw source code of the poisson-mfem.py example.
There is also a [commented version](@ref poisson-mfem-py)
of this example in the @ref examples section.

@include{lineno} poisson-mfem.py
