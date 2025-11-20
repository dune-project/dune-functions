@page poisson-mfem-py Mixed FE discretization of the Poisson equation (Python)
<!--
SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file AUTHORS.md
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception or LGPL-3.0-or-later
-->

The following example explains how to solve a mixed formulation of the Poisson problem
with Lagrange elements for the solution and Raviart-Thomas elements for the gradient.
The full example file @ref source-poisson-mfem-py can be found in the `examples/` subdirectory.

Be warned that writing a quadrature loop for a stiffness matrix with [dune-functions][]
and Python is currently quite slow. This approach is therefore mainly useful
for small problems, for educational purposes, or if you need complete control.
For faster and more convenient ways to assemble finite element stiffness matrices,
have a look at the [dune-fufem](https://gitlab.dune-project.org/fufem/dune-fufem)
and [dune-fem](https://gitlab.dune-project.org/dune-fem/dune-fem) modules.

## Mixed formulation of the Poisson problem

Let \f$ \Omega \f$ be a bounded open subset of \f$ \mathbb{R}^d \f$.
The Poisson problem asks for a scalar function \f$ u : \Omega \to \mathbb{R} \f$
that solves
\f[
 -\Delta u = f
 \qquad
 \text{on $\Omega$}
\f]
for a given function \f$f : \Omega \to \mathbb{R} \f$. In this example, we consider
boundary conditions
\f{alignat*}{{2}
 u & = u_0
 \qquad
 & \text{on $\Gamma_D$}
 \\
 \langle \nabla u, n \rangle & = g
 & \text{on $\Gamma_N$},
\f}
where \f$ \Gamma_D \f$ and \f$ \Gamma_N \f$ form a partition of the
domain boundary \f$ \partial \Omega \f$.

The mixed formulation of this problem introduces the gradient of \f$ u \f$
as a new variable \f$ \sigma \f$.  One obtains a system of two first-order equations
for two unknowns
\f{align*}
 \sigma - \nabla u & = 0
 \\
 -\operatorname{div} \sigma & = f,
\f}
and the Neumann boundary condition turns into
\f[
 \langle \sigma, n \rangle = g
 \qquad
 \text{on $\Gamma_N$.}
\f]
Multiplying with test functions \f$ \tau \f$ and \f$ v \f$ (where the normal component
of \f$ \tau \f$ is supposed to be zero \f$ \Gamma_n \f$) yields
\f{align*}
 \int_\Omega \langle \sigma, \tau \rangle\,dx
  - \int_\Omega \langle \nabla u, \tau \rangle \,dx & = 0 \\
 \int_\Omega \operatorname{div} \sigma \cdot v\,dx &= -\int_\Omega fv\,dx.
\f}

To obtain a symmetric bilinear form one applies integration by parts
to the first equation. One gets:
Find \f$ \sigma \in H_D(\text{div}) \f$ and \f$ u \in L^2(\Omega) \f$ such that
\f{alignat*}{{2}
 \int_\Omega \langle \sigma, \tau \rangle\,dx + \int_\Omega \langle \operatorname{div} \tau, u \rangle\,dx
 &=
 \int_{\Gamma_D} \langle \tau, n \rangle u\,ds
 & \qquad & \forall \tau \in H_0(\text{div})\\
 \int_\Omega \operatorname{div} \sigma \cdot v\,dx \qquad \qquad  & = - \int_\Omega fv\,dx && \forall v \in L^2(\Omega).
\f}
Here, the space \f$ H(\text{div}) \f$ is defined as
\f[
 H(\text{div}) := \Big\{ \tau \in L^2(\Omega)^d \mid \operatorname{div} \tau \in L^2(\Omega) \Big\}.
\f]
It is a vector space of vector-valued functions, and a strict
superset of \f$ H^1(\Omega)^d \f$.
The space \f$ H_D(\text{div}) \f$ denotes the subspace of \f$ H(\text{div}) \f$
of vector fields complying with the boundary condition on \f$ \Gamma_N \f$.
The space \f$ H_D(\text{div}) \f$ denotes the subspace of \f$ H(\text{div}) \f$
of vector fields having zero normal trace on \f$ \Gamma_N \f$.


A natural approximation space for \f$ H(\text{div}) \f$ is the space of Raviart-Thomas elements.
It consists of particular piecewise polynomial vector fields whose
normal components are continuous across element boundaries, but whose
tangential components may jump.  The solution \f$ u \f$ lives in \f$ L^2 \f$,
and therefore piecewise constant finite elements can be used.
Even though the boundary condition for \f$ u \f$ is not explicitly enforced, it can be shown
to hold in the limit of vanishing element diameter.

## The implementation

The Python implementation starts by importing required external modules.
It needs code from `numpy` and `scipy` for the sparse linear algebra,
and from `dune` for grids, finite element bases, and quadrature rules.

@snippet poisson-mfem.py imports

### Main program

The main program first creates the grid.  For simplicity, the domain here
is the unit square in two space dimensions, and the grid is a 50 by 50 quadrilateral grid:

@snippet poisson-mfem.py createGrid

The next line creates the bases of the two finite element spaces:
The Raviart-Thomas basis for \f$ \sigma \f$, and a zero-order Lagrange basis (i.e., piecewise
constant functions) for \f$ u \f$. In [dune-functions], such pairs of bases
are represented in a tree data structure with one root and two leaves:
Each of the leaves represents one of the two bases, and the root combines them into a single object.
This is the purpose of the `functions.composite` method in the following
code snippet. Setting up the basis tree also allows to select the way
how global degrees of freedom are numbered. Here, we use the default settings,
which lead to plain integers (actually, arrays of length 1) being used.

@snippet poisson-mfem.py createBasis

The method `subspaceBasis` gives access to the leaves of
a composite basis. The numbers 0 and 1 passed as arguments single out
the desired leaves. Note that accessing the two bases as subspace bases
of a tree is not the same as simply having two separate basis objects:
The actual basis functions are not affected by this, but the numbers differ.

Then, the problem stiffness matrix and right-hand side vector are assembled.
The code uses two separate methods for that, which are discussed below.
The source term \f$ f \f$ is specified in a lambda expression `f`, which here
encodes the constant function \f$ f : x \mapsto 2 \f$.

@snippet poisson-mfem.py assembly

The next code block manages the handling of the Dirichlet boundary conditions.
In the dual formulation used here, Dirichlet values are only applied
to the gradient field \f$ \sigma \f$, and only to its normal components.
This is easy to implement with Raviart-Thomas elements, whose degrees
of freedom are exactly the normal components of the vector field.

In this example, the gradient of \f$ u \f$ is prescribed on the upper
and lower horizontal edge of the domain boundary.
To find all degrees of freedom on these grid boundary parts, the method
`markBoundaryDOFsByIndicator` fills the array `isDirichlet` with the
corresponding information: The array contains a boolean variable per
degree of freedom, which is set to `True` if a degree of freedom is attached
to a grid intersection whose center is on the relevant boundary.

@snippet poisson-mfem.py dirichletDOFs

Then, the function \f$ g \f$ from the Neumann boundary condition
\f$ \langle \sigma, n \rangle = g \f$ is defined. In this example,
it is given by the closed-form expression
\f[
 g(x_0,x_1) = \sin(2\pi x_0).
\f]
This expression is evaluated
at the boundary Lagrange nodes by the method `basis.interpolate`.

\todo Remove the need to involve the normal vector here!

@snippet poisson-mfem.py dirichletValues

After this, the information about the boundary condition is integrated into
the linear system of equations. This happens by calling a method that appears
earlier in the file:

@snippet poisson-mfem.py dirichletIntegration

The implementation of this method is

@snippet poisson-mfem.py dirichletIntegration

The code loops over all nonzero entries of the
matrix. If an entry belongs to a matrix row that corresponds to a boundary
degree of freedom, it is set to 1 if it is the diagonal entry, and to 0 otherwise.
This pins the degree of freedom to the corresponding boundary value.
This part of the code is completely independent from %Dune.

At this point, the linear system of equations has been completely assembled,
and can now be solved. For this, the code simply calls a sparse direct solver
from the Python `scipy` package.

@snippet poisson-mfem.py solving

The solver leaves the algebraic solution vector in the vector variable `x`.
The solution is then written to a file, in the VTK file format.  To do so,
the solution vector must be reinterpreted as a finite element functions
again. The `asFunction` method does this:

@snippet poisson-mfem.py vtkWriting

This will write a file called `poisson-mfem-result.vtu`, which contains the grid
and the solution \f$ u \f$ in a field called "value", and the approximate
gradient vector field \f$ \sigma \f$ in a field called "gradient". It can be opened, for example,
with the [ParaView](https://www.paraview.org/) visualization program.

## The assembler

The assembler is implemented in four methods above the main method of the program.
Two of them implement the local assembler, i.e., they compute stiffness matrix
and load vector for a single element. The two global assembler methods
then piece together the local matrices and vectors into a single global matrix
and vector, respectively.

### Global assembler

The function that assembles the global stiffness matrix is `assembleMixedPoissonMatrix`.
As described in many text books, this assembly is implemented as a loop over the grid
elements. On each element, all relevant integrals are computed and stored
in a matrix `elementMatrix`, the *element stiffness matrix*. Since the basis functions
of both finite element spaces are zero on most elements by construction, this element
stiffness matrix is much smaller than the global stiffness matrix,
but it is dense. Computation of `elementMatrix` happens by calling the method
`getLocalMatrix` (discussed below). Then, the entries of the element
stiffness matrix are added to the appropriate places of the global
`stiffnessMatrix` object. The `basis` object provides the correspondence
between the entries of the element stiffness matrices and the entries
of the global stiffness matrix. Note how this does not require to distinguish
between value and gradient degrees of freedom at all.

@snippet poisson-mfem.py assembleMixedPoissonMatrix

Note that `stiffnessMatrix` is sparse (i.e., it contains mainly zeros). This property
has to be exploited. Therefore, a `lil_matrix` object from the Python `scipy`
package is used to build it. After assembly has been completed, the matrix is
converted to Compressed Sparse Row (CSR) format, because that is more
efficient for the solver.

The next method does the global assembly for the load vector.
Since a nonzero volume term only appears for the second equation, only the
second basis (accessed via `functions.subspaceBasis(basis,1)`) is used.
The volume term function given by the main method takes as argument
a point in world coordinates. However, the local assembler needs a function
that takes arguments in element-local coordinates. The C++ interface
provides the class Dune::Functions::AnalyticGridViewFunction to do
this conversion, but this is not exposed as Python yet.
Instead, we project the `volumeTerm` onto the Lagrange finite element space.
This serves the same purpose.

@snippet poisson-mfem.py assembleMixedPoissonRhs

### Local assembler

The element stiffness matrices are computed by a free function called `getLocalMatrix`.
Let \f$ (\theta_i) \f$ be the Raviart-Thomas shape functions and \f$ (\varphi_j) \f$
the zero-order Lagrange shape functions. The element and the set of basis functions on the element
are passed to the method via a `localView` object. See the Dune::Functions::DefaultLocalView page
for documentation of the corresponding C++ interface.

For each element \f$ T \f$, the stiffness matrix has a two-by-two block form
\f[
 (A_T) = \begin{pmatrix} C & B \\ B^T & 0 \end{pmatrix},
\f]
with entries
\f[
 C_{ij} = \int_\Omega \langle \theta_i, \theta_j\,dx
 \qquad
 B_{ij} = \int_\Omega \operatorname{div} \theta_i, \phi_j \rangle\,dx
\f]
Both submatrices are dense objects. In the code, they are referred to as
"gradient-gradient coupling" and "value-gradient coupling", respectively.

The computation of the integrals happens via numerical quadrature, with a
quadrature rule that is defined on a reference element \f$ T_\text{ref} \f$.
The transformation of the divergence uses the transformation of the full
Jacobian \f$ \nabla \theta_i \f$ as explained in the @ref poisson-pq2-py example.

@snippet poisson-mfem.py getLocalMatrix

The final method computes the load vector for one element. Since only
the second equation has a nontrivial linear term, the entries corresponding
to Raviart-Thomas degrees of freedom are zero. The entries for the Lagrange
degrees of freedom are
\f[
 (b_T)_i = \int_T f \varphi_i\,dx.
\f]

The element load vector is computed by the following code:

@snippet poisson-mfem.py getLocalVolumeTerm


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
