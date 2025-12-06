# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#! [imports]
import numpy as np
import scipy.sparse.linalg
from scipy.sparse import lil_matrix

import dune.geometry
import dune.grid as grid
import dune.functions as functions
#! [imports]

# Compute element stiffness matrix and element load vector
#
# TODO: This assembler loop is very inefficient in terms of run time and should be improved using Python vectorization.
# See discussion at https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/295 for hints and code pointers.
#! [localAssembler]
def localAssembler(localView, volumeTerm):

    # Number of degrees of freedom on this element
    n = len(localView)

    # The grid element
    element = localView.element()

    # The set of shape functions
    localBasis = localView.tree().finiteElement().localBasis

    # Initialize element matrix and load vector
    localA = np.zeros((n,n))
    localB = np.zeros(n)

    # Create a quadrature rule and integrate
    quadRule = dune.geometry.quadratureRule(element.type, order=4)
    for pt in quadRule:

        # The determinant that appears in the integral transformation formula
        integrationElement = element.geometry.integrationElement(pt.position)

        # Evaluate all shape functions (The return value is an array!)
        values = localBasis.evaluateFunction(pt.position)

        # Evaluate the shape function Jacobians on the reference element (array of arrays)
        referenceJacobians = localBasis.evaluateJacobian(pt.position)

        # Transform the reference Jacobians to the actual element
        geometryJacobianInverse = element.geometry.jacobianInverse(pt.position)
        jacobians = [ np.dot(np.array(g)[0], geometryJacobianInverse) for g in referenceJacobians ]

        quadPosGlobal = element.geometry.toGlobal(pt.position)

        for i in range( n ):
            for j in range( n ):
                localA[i,j] += pt.weight * integrationElement * np.dot(jacobians[i], jacobians[j])

            localB[i] += pt.weight * integrationElement * values[i] * volumeTerm(quadPosGlobal)

    return localA, localB
#! [localAssembler]


# The assembler for the global stiffness matrix
#! [assembleLaplaceMatrix]
def assembleLaplaceMatrix(basis, volumeTerm):

    # Total number of degrees of freedom
    n = len(basis)

    # Make empty sparse matrix
    A = lil_matrix((n,n))

    # Make empty vector
    b = np.zeros(n)

    # View on the finite element basis on a single element
    localView = basis.localView()

    # Loop over all grid elements
    grid = basis.gridView
    for element in grid.elements:

        # Bind the localView to the current element
        localView.bind(element)

        # Number of degrees of freedom on the current element
        localN = len(localView)

        # Assemble the local stiffness matrix and load vector
        localA, localb = localAssembler(localView, volumeTerm)

        # Copy the local entries into the global matrix and vector
        for i in range(localN):

            gi = localView.index(i)[0]

            b[gi] += localb[i]

            for j in range(localN):
                gj = localView.index(j)[0]
                A[gi, gj] += localA[i, j]

    # Convert matrix to CSR format
    return A.tocsr(), b
#! [assembleLaplaceMatrix]


############################  main program  ###################################

#! [createGrid]
# Number of grid elements in one direction
gridSize = 32

# Create a grid of the unit square
grid = grid.structuredGrid([0,0],[1,1],[gridSize,gridSize])
#! [createGrid]

#! [createBasis]
# Create a second-order Lagrange FE basis
basis = functions.defaultGlobalBasis(grid, functions.Lagrange(order=2))
#! [createBasis]

#! [assembly]
# Source term
f = lambda x : 10

# Compute stiffness matrix and load vector
A,b = assembleLaplaceMatrix(basis, f)
#! [assembly]

#! [dirichletDOFs]
# Determine all coefficients that are on the boundary
isDirichlet = np.zeros(len(basis))

def markDOF(boundaryDOFNumber):
    isDirichlet[boundaryDOFNumber] = True

functions.boundarydofs.forEachBoundaryDOF(basis,markDOF)
#! [dirichletDOFs]

#! [dirichletValues]
# The function that implements the Dirichlet values
dirichletValueFunction = lambda x : np.sin(2*np.pi*x[0])

# Get coefficients of a Lagrange-FE approximation of the Dirichlet values
dirichletValues = np.zeros(len(basis))
basis.interpolate(dirichletValues, dirichletValueFunction)
#! [dirichletValues]

#! [dirichletIntegration]
# Integrate Dirichlet conditions into the matrix and load vector
rows, cols = A.nonzero()

for i,j in zip(rows, cols):
    if isDirichlet[i]:
        if i==j:
            A[i,j] = 1.0
        else:
            A[i,j] = 0
        b[i] = dirichletValues[i]
#! [dirichletIntegration]

#! [solving]
# Solve linear system!
x = scipy.sparse.linalg.spsolve(A, b)
#! [solving]

#! [vtkWriting]
# Write result as vtu file
uh = basis.asFunction(x)

vtkWriter = grid.vtkWriter(subsampling=2)
uh.addToVTKWriter("u", vtkWriter, dune.grid.DataType.PointData)
vtkWriter.write("poisson-pq2-result")
#! [vtkWriting]
