@page poisson-pq2-cpp Poisson equation
<!--
SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file AUTHORS.md
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception or LGPL-3.0-or-later
-->

The following explains how to solve the Poisson equation
using [dune-functions][]. The full example @ref source-poisson-pq2-cpp
can be found in the `examples/` subdirectory.



### Local assemblers

The program first includes the required headers and then defines
free functions for assembling the Poisson problem.
The function `getLocalMatrix()` implements the assembler
of the local stiffness matrix for the bilinear form
@f$(u,v) \mapsto \int_\Omega \nabla u(x)\nabla v(x)dx@f$.

@snippet poisson-pq2.cc getLocalMatrix

The `getVolumeTerm()` functions implements the local assembler
for the volume right hand side term @f$\int_\Omega f(x)v(x)dx@f$.

@snippet poisson-pq2.cc getVolumeTerm


### Global assembler

Assemble the global matrix pattern.

@snippet poisson-pq2.cc getOccupationPattern


Assembly of matrix and right-hand-side.

@snippet poisson-pq2.cc assembleLaplaceMatrix


### Helper functions

Treatment of boundary condition.

@snippet poisson-pq2.cc boundaryTreatment

Create a mixed grid containing triangles and quadrilaterals.

@snippet poisson-pq2.cc createMixedGrid



### Setup and initialization

Initialize MPI

@snippet poisson-pq2.cc startMainAndMPI



### Create grid and finite element basis

Create a mixed grid and obtain a leaf grid view.

@note
A `GridView` is a view to a subset of the grid's elements, vertices, ...
that should be stored by value and can be copied cheaply.
Grids in [dune][] are in general hierarchical and composed by elements
on several levels. The discretization usually lives on the set of
most refined elements that is denote the leaf `GridView` in dune.

@snippet poisson-pq2.cc createGridCall

As a next step the program creates a global finite element
basis on the `GridView`.

@snippet poisson-pq2.cc makeBasis

And here is the rest of the file:

@snippet poisson-pq2.cc theRest


[dune]: https://dune-project.org
[dune-functions]: https://gitlab.dune-project.org/staging/dune-functions


<div class="section_buttons">
| Previous      |                         Next |
|:--------------|-----------------------------:|
| @ref examples |                              |
</div>



@page source-poisson-pq2-cpp poisson-pq2.cc

This is the raw source code of the poisson-pq2.cc example.
There is also a [commented version](@ref poisson-pq2-cpp)
of this example in the @ref examples section.

@include{lineno} poisson-pq2.cc
