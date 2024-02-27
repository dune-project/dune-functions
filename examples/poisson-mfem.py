import numpy
import scipy.sparse.linalg
from scipy.sparse import lil_matrix
from io import StringIO

import dune.geometry
import dune.grid as grid
import dune.functions as functions

# Compute the stiffness matrix for a single element
def getLocalMatrix(localView):

    n = len(localView)    # Number of local degrees of freedom (flux + pressure)

    # Get the grid element from the local FE basis view
    element = localView.element()

    # Make dense element stiffness matrix
    elementMatrix = numpy.zeros((n,n))

    # Get set of shape functions for this element
    fluxLocalFiniteElement     = localView.tree().child(0).finiteElement()
    pressureLocalFiniteElement = localView.tree().child(1).finiteElement()

    # The actual shape functions on the reference element
    localFluxBasis = fluxLocalFiniteElement.localBasis
    localPressureBasis = pressureLocalFiniteElement.localBasis

    nFlux = len(fluxLocalFiniteElement)
    nPressure = len(pressureLocalFiniteElement)

    # Get a quadrature rule
    fluxOrder = dim*fluxLocalFiniteElement.localBasis.order
    pressureOrder = dim*pressureLocalFiniteElement.localBasis.order
    quadOrder = numpy.max((2*fluxOrder, (fluxOrder-1)+pressureOrder))

    quadRule = dune.geometry.quadratureRule(element.type, quadOrder)

    # Loop over all quadrature points
    for quadPoint in quadRule:

        # Position of the current quadrature point in the reference element
        quadPos = quadPoint.position

        # The inverse Jacobian of the map from the reference element to the element
        geometryJacobianInverse = element.geometry.jacobianInverse(quadPos)

        # The multiplicative factor in the integral transformation formula
        integrationElement = element.geometry.integrationElement(quadPos)

        # --------------------------------------------------------------
        #  Shape functions - flux
        # --------------------------------------------------------------

        # Values of the flux shape functions on the current element
        fluxValues = localFluxBasis.evaluateFunction(quadPos)

        # Gradients of the flux shape function gradients on the reference element
        fluxReferenceJacobians = localFluxBasis.evaluateJacobian(quadPos)

        fluxDivergence = numpy.zeros(nFlux)

        # Domain transformation of Jacobians and computation of div = trace(Jacobian)
        # TODO: Extend the Dune Python interface to allow to do this without casting to numpy.array
        for i in range(nFlux):
            fluxDivergence[i] = numpy.trace(numpy.array(fluxReferenceJacobians[i]) * geometryJacobianInverse)

        # --------------------------------------------------------------
        #  Shape functions - pressure
        # --------------------------------------------------------------

        # Values of the pressure shape functions
        pressureValues = localPressureBasis.evaluateFunction(quadPos)

        # --------------------------------------------------------------
        #  Flux--flux coupling
        # --------------------------------------------------------------

        for i in range(nFlux):

            row = localView.tree().child(0).localIndex(i)

            for j in range(nFlux):

                col = localView.tree().child(0).localIndex(j)
                elementMatrix[row,col] += numpy.dot(fluxValues[i], fluxValues[j]) * quadPoint.weight * integrationElement

        # --------------------------------------------------------------
        #  Flux--pressure coupling
        # --------------------------------------------------------------

        for i in range(nFlux):

            fluxIndex = localView.tree().child(0).localIndex(i)

            for j in range(nPressure):

                pressureIndex = localView.tree().child(1).localIndex(j)

                # Pre-compute matrix contribution
                tmp = (fluxDivergence[i] * pressureValues[j][0]) * quadPoint.weight * integrationElement

                elementMatrix[fluxIndex, pressureIndex] += tmp
                elementMatrix[pressureIndex, fluxIndex] += tmp

    return elementMatrix


# Compute the right-hand side for a single element
def getVolumeTerm(localView, localVolumeTerm):

    # Get the grid element from the local FE basis view
    element = localView.element()

    n = len(localView)
    localRhs = numpy.zeros(n)

    # Get set of shape functions for this element
    # Only the pressure part has a non-zero right-hand side
    pressureLocalFiniteElement = localView.tree().child(1).finiteElement()

    # A quadrature rule
    dim = element.dimension
    quadOrder = 2*dim*pressureLocalFiniteElement.localBasis.order
    quadRule = dune.geometry.quadratureRule(element.type, quadOrder)

    nPressure = len(pressureLocalFiniteElement)

    # Loop over all quadrature points
    for quadPoint in quadRule:

        # Position of the current quadrature point in the reference element
        quadPos = quadPoint.position

        # The multiplicative factor in the integral transformation formula
        integrationElement = element.geometry.integrationElement(quadPos)

        # Evaluate the strong right-hand side at the quadrature point
        functionValue = localVolumeTerm(quadPos)

        # Evaluate all shape function values at this point
        pressureValues = pressureLocalFiniteElement.localBasis.evaluateFunction(quadPos)

        # Actually compute the vector entries
        for j in range(nPressure):
            pressureIndex = localView.tree().child(1).localIndex(j)
            localRhs[pressureIndex] += - pressureValues[j][0] * functionValue * quadPoint.weight * integrationElement

    return localRhs


# Assemble the divergence stiffness matrix on the given grid view
def assembleMixedPoissonMatrix(basis):

    # Get the grid view from the finite element basis
    gridView = basis.gridView

    n = len(basis)

    # Make an empty stiffness matrix
    stiffnessMatrix = lil_matrix( (n,n) )

    # A view on the FE basis on a single element
    localView = basis.localView()

    # A loop over all elements of the grid
    for element in basis.gridView.elements:

        # Bind the local FE basis view to the current element
        localView.bind(element)

        # Now let's get the element stiffness matrix
        # A dense matrix is used for the element stiffness matrix
        elementMatrix = getLocalMatrix(localView)

        # Add element stiffness matrix onto the global stiffness matrix
        for i in range(len(localView)):

            # The global index of the i-th local degree of freedom of the element
            row = localView.index(i)[0]

            for j in range(len(localView)):

                # The global index of the j-th local degree of freedom of the element
                col = localView.index(j)[0]
                stiffnessMatrix[row,col] += elementMatrix[i, j];

    # Transform the stiffness matrix to CSR format, and return it
    return stiffnessMatrix.tocsr()


# Assemble the divergence stiffness matrix on the given grid view
def assembleMixedPoissonRhs(basis, volumeTerm):

    # Get the grid view from the finite element basis
    gridView = basis.gridView

    # Get the basis for the pressure variable
    pressureBasis = functions.subspaceBasis(basis, 1)

    # Represent the volume term function as a FE function in the pressure space
    volumeTermCoeff = numpy.zeros(len(basis))
    pressureBasis.interpolate(volumeTermCoeff, volumeTerm)
    volumeTermGF = pressureBasis.asFunction(volumeTermCoeff)

    # A view on a single element
    localVolumeTerm = volumeTermGF.localFunction()

    # Set rhs to correct length -- the total number of basis vectors in the basis
    n = len(basis)
    rhs = numpy.zeros(n)

    # A view on the FE basis on a single element
    localView = basis.localView()

    # A loop over all elements of the grid
    for element in gridView.elements:

        # Bind the local FE basis view to the current element
        localView.bind(element)

        # Now get the local contribution to the right-hand side vector
        localVolumeTerm.bind(element)
        localRhs = getVolumeTerm(localView, localVolumeTerm)

        for i in range(len(localRhs)):
            # The global index of the i-th vertex of the element
            row = localView.index(i)[0]
            rhs[row] += localRhs[i]

    return rhs


# Mark all DOFs associated to entities for which # the boundary intersections center
# is marked by the given indicator functions.
#
# This method simply calls the corresponding C++ code.  A more Pythonic solution
# is planned to appear eventually...
def markBoundaryDOFsByIndicator(basis, vector, indicator):
    code="""
    #include<utility>
    #include<functional>
    #include<dune/common/fvector.hh>
    #include<dune/functions/functionspacebases/boundarydofs.hh>
    template<class Basis, class Vector, class Indicator>
    void run(const Basis& basis, Vector& vector, const Indicator& indicator)
    {
      auto vectorBackend = vector.mutable_unchecked();
      Dune::Functions::forEachBoundaryDOF(basis, [&] (auto&& localIndex, const auto& localView, const auto& intersection) {
        if (indicator(intersection.geometry().center()).template cast<bool>())
          vectorBackend[localView.index(localIndex)] = true;
      });
    }
    """
    dune.generator.algorithm.run("run",StringIO(code), basis, vector, indicator)


# This incorporates essential constraints into matrix # and rhs of a linear system.
# The mask vector isConstrained # indicates which DOFs should be constrained,
# x contains the desired values of these DOFs. Other entries of x # are ignored.
# Note that this implements the symmetrized approach to modify the matrix.
def incorporateEssentialConstraints(A, b, isConstrained, x):
    b -= A*(x*isConstrained)
    N = len(b)
    rows, cols = A.nonzero()
    for i,j in zip(rows, cols):
        if isConstrained[i] or isConstrained[j]:
          A[i,j] = 0
    for i in range(N):
        if isConstrained[i]:
            A[i,i] = 1
            b[i] = x[i]



############################  main program  ###################################

# Number of grid elements per direction
dim = 2
elements = [50, 50]
l = [1, 1]

# Create a grid of the unit square
gridView = grid.structuredGrid([0,0],l,elements)

# Construct a pair of finite element space bases
# Note: In contrast to the corresponding C++ example we are using a single matrix with scalar entries,
# and plain numbers to index it (no multi-digit multi-indices).
k = 0  # order
basis = functions.defaultGlobalBasis(gridView, functions.Composite(functions.RaviartThomas(order=k),
                                                                             functions.Lagrange(order=k),
                                                                             blocked=False,
                                                                             layout="lexicographic"))

fluxBasis = functions.subspaceBasis(basis, 0);

pressureBasis = functions.subspaceBasis(basis, 1);

# Compute the stiffness matrix and the load vector
stiffnessMatrix = assembleMixedPoissonMatrix(basis)

# The volume source term
rightHandSide = lambda x : 2

rhs = assembleMixedPoissonRhs(basis, rightHandSide)

# This marks the top and bottom boundary of the domain
fluxDirichletIndicator = lambda x : 1.*((x[dim-1] > l[dim-1] - 1e-8) or (x[dim-1] < 1e-8))

############################################################
# ToDo: We should provide binding for FaceNormalGridFunction
# and support for grid functions to the bindings of interpolate().
# This would allow to avoid having to define the normal field
# manually.
############################################################

normal = lambda x : numpy.array([0.,1.]) if numpy.abs(x[1]-1) < 1e-8 else numpy.array([0., -1.])
fluxDirichletValues = lambda x : numpy.sin(2.*numpy.pi*x[0]) * normal(x)

isDirichlet = numpy.zeros(len(basis))

# Mark all DOFs located in a boundary intersection marked
# by the fluxDirichletIndicator function. If the flux
# ansatz space also contains tangential components, this
# approach will fail, because those are also marked.
# For Raviart-Thomas this does not happen.
markBoundaryDOFsByIndicator(fluxBasis, isDirichlet, fluxDirichletIndicator);

# ToDo: This should be constrained to boundary DOFs
fluxDirichletCoeffs = numpy.zeros(len(basis))
fluxBasis.interpolate(fluxDirichletCoeffs, fluxDirichletValues);

# //////////////////////////////////////////
#   Modify Dirichlet rows
# //////////////////////////////////////////

incorporateEssentialConstraints(stiffnessMatrix, rhs, isDirichlet, fluxDirichletCoeffs)

# //////////////////////////
#    Compute solution
# //////////////////////////

x = scipy.sparse.linalg.spsolve(stiffnessMatrix, rhs)

# ////////////////////////////////////////////////////////////////////////////////////////////
#   Write result to VTK file
# ////////////////////////////////////////////////////////////////////////////////////////////

# TODO: Improve file writing.  Currently this simply projects everything
# onto a first-order Lagrange space
vtkWriter = gridView.vtkWriter(2)

fluxFunction = fluxBasis.asFunction(x)
fluxFunction.addToVTKWriter("flux", vtkWriter, grid.DataType.PointVector)

pressureFunction = pressureBasis.asFunction(x)
pressureFunction.addToVTKWriter("pressure", vtkWriter, grid.DataType.PointData)

vtkWriter.write("poisson-mfem-result")
