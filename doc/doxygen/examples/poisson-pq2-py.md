@page poisson-pq2-py Poisson equation (Python)
<!--
SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file AUTHORS.md
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception or LGPL-3.0-or-later
-->

The following example explains how to solve the Poisson problem with Lagrange
finite elements using the Python interface of [dune-functions][].
The full example file @ref source-poisson-pq2-py can be found in the `examples/` subdirectory.

Be warned that writing a quadrature loop for a stiffness matrix with [dune-functions][]
and Python is currently quite slow. This approach is therefore mainly useful
for small problems, for educational purposes, or if you need complete control.
For faster and more efficient ways to assemble finite element stiffness matrices,
have a look at the [dune-fufem](https://gitlab.dune-project.org/fufem/dune-fufem)
and [dune-fem](https://gitlab.dune-project.org/dune-fem/dune-fem) modules.

## The Poisson problem

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
The finite element method replaces the partial differential equation by a so-called
*weak formulation*, and solves it in a space \f$ V_h \f$ of continuous,
piecewise polynomial functions: The goal is now to find a function \f$ u_h \in V_h \f$
still satisfying the boundary conditions, such that
\f[
 \int_\Omega \langle \nabla u_h, \nabla v_h \rangle \,dx
 =
 \int_\Omega f v_h\,dx
 \qquad
 \forall v_h \in V_{h,0},
\f]
where \f$ V_{h,0} \f$ is the space of finite element functions that are zero
on the boundary.
Writing this weak form in terms of a basis \f$ \{ \varphi_i \} \f$ of \f$ V_h \f$, one obtains a sparse
linear system of equations
\f[
 Ax = b,
\f]
which can be solved using standard algorithms from numerical linear algebra.

## The implementation

As every Python program, this one starts by importing required external modules.
In this case, we need code from `numpy` and `scipy` for the sparse linear algebra,
and from `dune` for grids, finite element bases, and quadrature rules.

@snippet poisson-pq2.py imports


### Main program

The main program starts by creating a grid.  For simplicity, the domain here
is the unit square in two space dimensions, and the grid is a 32 by 32 quadrilateral grid:

@snippet poisson-pq2.py createGrid

The next line creates the basis of the finite element space.  We use second-order
Lagrange finite elements here, i.e., the space \f$ V_h \f$ consists of continuous
functions that are quadratic in each variable on each element. The class
`functions.Lagrange` encodes the nodal basis of that space:

@snippet poisson-pq2.py createBasis

Then, the Poisson problem is assembled. The source term \f$ f \f$ is specified
in a lambda expression `f`, which here encodes the constant function \f$ f : x \mapsto 10 \f$.
The actual assembly then happens in the method `assembleLaplaceMatrix`, which we
discuss below.

@snippet poisson-pq2.py assembly

The `assembleLaplaceMatrix` method returns *two* objects, the stiffness matrix `A`
and the assembled source term vector `b`.

The next code block manages the handling of the Dirichlet boundary conditions.
In this example, the entire boundary is the Dirichet boundary.
To find all degrees of freedom on the grid boundary, the method `forEachBoundaryDOF`
visits each degree of freedom (DOF) associated to the grid boundary and
calls the callback function `markDOF` with the global index of that degree of freedom.
This callback then sets the corresponding entry in an array
called `isDirichlet` to `True`.

@snippet poisson-pq2.py dirichletDOFs

Then, the function \f$ g \f$ from the problem formulation
is defined. In this example the Dirichlet values are given by the closed-form
expression
\f[
 g(x_0,x_1) = \sin(2\pi x_0).
\f]
This closed form expression is evaluated
at the boundary Lagrange nodes by the method `basis.interpolate`.

@snippet poisson-pq2.py dirichletValues

After this, the information about the boundary conditions is integrated into
the linear system of equations. The code loops over all nonzero entries of the
matrix. If an entry belongs to a matrix row that corresponds to a boundary
degree of freedom, it is set to 1 if it is the diagonal entry, and to 0 otherwise.
This pins the degree of freedom to the corresponding boundary value.
This part of the code is completely independent from %Dune.

@snippet poisson-pq2.py dirichletIntegration

At this point, the linear system of equations has been completely assembled,
and can now be solved. For this, the code simply calls a sparse direct solver
from the Python `scipy` package. This is not part of %Dune, either:

@snippet poisson-pq2.py solving

The solver leaves the algebraic solution vector in the vector variable `x`.
The solution is then written to a file, in the VTK file format.  To do so,
the solution vector must be reinterpreted as a finite element functions
again. The `asFunction` method does this:

@snippet poisson-pq2.py vtkWriting

This will write a file called `poisson-pq2-result.vtu`, which contains the grid
and the solution function \f$ u_h \f$. It can be opened, for example,
with the [ParaView](https://www.paraview.org/) visualization program.

## The assembler

The assembler is implemented in two methods above the main method of the program.

### Global assembler
The central function is `assembleLaplaceMatrix`, which computes the stiffness matrix \f$ A \f$
with entries
\f[
 A_{ij} = \int_\Omega \langle \nabla \varphi_i, \nabla \varphi_j \rangle dx
\f]
and the load vector \f$ b \f$ with entries
\f[
 b_i = \int_\Omega f \varphi_i\,dx.
\f]
As described in many text books, this assembly is implemented as a loop over the grid
elements. On each element, all relevant integrals are computed and stored
in a matrix `localA`, the *element stiffness matrix*. Since the basis functions
\f$ \varphi_i \f$ are zero on most elements by construction, this element
stiffness matrix is much smaller than the global stiffness matrix `A`,
but it is dense. Computation of `localA` happens by calling yet another method
`localAssembler` (discussed below). Afterwards, the entries of the element
stiffness matrix are added to the appropriate places of the global
stiffness matrix `A`. The `basis` object provides the correspondence
between the entries of the element stiffness matrices and the entries
of the global stiffness matrix.

@snippet poisson-pq2.py assembleLaplaceMatrix

Note that `A` is sparse (i.e., it contains mainly zeros). This property
has to be exploited. Therefore, a `lil_matrix` object from the Python `scipy`
package is used to build it. After assembly has been completed, the matrix is
converted to Compressed Sparse Row (CSR) format, because that is more
efficient for the solver.

### Local assembler

The element stiffness matrices are computed by a free function called `localAssembler`.
For a given element \f$ T \f$, it computes the element stiffness matrix
\f[
 (A_T)_{ij} = \int_T \langle \nabla \varphi_i, \nabla \varphi_j \rangle \,dx
\f]
and the element load vector \f$ b_T \f$ with entries
\f[
 (b_T)_i = \int_T f \varphi_i\,dx.
\f]
Both are dense objects. In this example, where the domain is two-dimensional,
the grid consists of quadrilaterals, and the finite element space consists
of second-order Lagrange elements, each \f$ A_T \f$ is a \f$ 9 \times 9\f$
matrix. The element and the set of basis functions on the element
are passed to the method via a `localView` object. See the @ref Dune::Functions::DefaultLocalView page
for documentation of the corresponding C++ interface.

The computation of the integrals happens via numerical quadrature, and the
quadrature rule is defined on a reference element \f$ T_\text{ref} \f$.
The integral therefore needs to be transformed to the reference element,
using a mapping \f$ F : T_\text{ref} \to T \f$:
\f[
 (A_T)_{pq}
  =
  \int_{T_\text{ref}}
  \Big \langle (\nabla_\xi \hat{\varphi}_{T,p}(\xi)) \cdot (\nabla_\xi F_T(\xi))^{-1},
  (\nabla_\xi \hat{\varphi}_{T,q}(\xi)) \cdot (\nabla_\xi F_T(\xi))^{-1} \Big\rangle
  | \det \nabla_\xi F_T(\xi)|\,d\xi.
\f]
The integrals over \f$ T_\text{ref} \f$ are then computed by a quadrature rule.


@snippet poisson-pq2.py localAssembler


[dune]: https://dune-project.org
[dune-functions]: https://gitlab.dune-project.org/staging/dune-functions


<div class="section_buttons">
| Previous      |                         Next |
|:--------------|-----------------------------:|
| @ref poisson-pq2-cpp | @ref poisson-mfem-py  |
</div>



@page source-poisson-pq2-py poisson-pq2.py

This is the raw source code of the poisson-pq2.py example.
There is also a [commented version](@ref poisson-pq2-py)
of this example in the @ref examples section.

@include{lineno} poisson-pq2.py
