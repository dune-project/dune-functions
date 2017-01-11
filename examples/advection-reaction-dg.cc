// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>
#include <cmath>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

/** \file
 * \brief Example implementation of a DG-discretized advection-reaction equation
 *
 * Implemented using Alexandre Ern's lecture notes "Discontinuous Galerkin Methods", Chapter 3.3.
 */

using namespace Dune;

// Compute the element stiffness matrix for an element coupling with itself
template <class LocalView, class MatrixType, class LocalVelocityField, class LocalReactionCoefficient>
void getLocalMatrix(const LocalView& localView,
                    MatrixType& elementMatrix,
                    LocalVelocityField&& localVelocityField,
                    LocalReactionCoefficient&& localReactionCoefficient)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();

  // Set all matrix entries to zero
  elementMatrix.setSize(localFiniteElement.size(),localFiniteElement.size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // Get a quadrature rule
  int order = 2*dim*localFiniteElement.localBasis().order();
  const auto& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad)
  {
    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPoint.position());

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPoint.position());

    // The values of the shape functions at the quadrature point
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(quadPoint.position(), values);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPoint.position(), referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
    for (size_t i=0; i<gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);

    // Compute the actual matrix entries
    for (size_t i=0; i<elementMatrix.N(); i++)
      for (size_t j=0; j<elementMatrix.M(); j++ )
      {
        // First: the reaction part
        elementMatrix[i][j] += localReactionCoefficient(quadPoint.position()) * values[i] * values[j] * quadPoint.weight() * integrationElement;

        elementMatrix[i][j] += ( localVelocityField(quadPoint.position()) * gradients[i] ) * values[j] * quadPoint.weight() * integrationElement;
      }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  //  Assemble the skeleton parts
  ////////////////////////////////////////////////////////////////////////////////////////////

  for (auto&& intersection : intersections(localView.globalBasis().gridView(), element))
  {
    // Get a quadrature rule on the element face
    int order = (dim-1)*localFiniteElement.localBasis().order();
    const auto& quad = QuadratureRules<double, dim-1>::rule(intersection.type(), order);

    for (auto&& quadPoint : quad)
    {
      auto positionInElement = intersection.geometryInInside().global(quadPoint.position());

      // The values of the shape functions at the quadrature point
      std::vector<FieldVector<double,1> > values;
      localFiniteElement.localBasis().evaluateFunction(positionInElement, values);

      auto factor = localVelocityField(positionInElement)*intersection.integrationOuterNormal(quadPoint.position());

      // Disable the term at outflow boundary faces
      if (!intersection.neighbor())
        factor = std::min(factor, 0.0);

      for (size_t i=0; i<values.size(); i++)
        for (size_t j=0; j<values.size(); j++)
          elementMatrix[i][j] += -1 * factor * values[i] * values[j] * quadPoint.weight();
    }

  }

}

// Compute the element stiffness matrix for an element coupling with a neighbor
template <class Intersection, class LocalView, class MatrixType, class LocalVelocityField>
void getOffDiagonalLocalMatrix(const Intersection& intersection,
                               const LocalView& insideLocalView,
                               const LocalView& outsideLocalView,
                               MatrixType& elementMatrix,
                               LocalVelocityField&& localVelocityField)
{
  const int intersectionDim = Intersection::mydimension;

  // Get set of shape functions for the inside and the outside element
  const auto& insideLocalFiniteElement = insideLocalView.tree().finiteElement();
  const auto& outsideLocalFiniteElement = outsideLocalView.tree().finiteElement();

  // Set all matrix entries to zero
  elementMatrix.setSize(insideLocalFiniteElement.size(),outsideLocalFiniteElement.size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // Get a quadrature rule on the element face
  int order = insideLocalFiniteElement.localBasis().order() * outsideLocalFiniteElement.localBasis().order();
  const auto& quad = QuadratureRules<double, intersectionDim>::rule(intersection.type(), order);

  for (auto&& quadPoint : quad)
  {
    auto positionInInsideElement = intersection.geometryInInside().global(quadPoint.position());
    auto positionInOutsideElement = intersection.geometryInOutside().global(quadPoint.position());

    // The values of the shape functions at the quadrature point
    std::vector<FieldVector<double,1> > insideValues, outsideValues;
    insideLocalFiniteElement. localBasis().evaluateFunction(positionInInsideElement,  insideValues);
    outsideLocalFiniteElement.localBasis().evaluateFunction(positionInOutsideElement, outsideValues);

    // Velocity field times face normal
    auto factor = localVelocityField(positionInInsideElement)*intersection.integrationOuterNormal(quadPoint.position());

    for (size_t i=0; i<insideValues.size(); i++)
      for (size_t j=0; j<outsideValues.size(); j++)
        elementMatrix[i][j] += factor * insideValues[i] * outsideValues[j] * quadPoint.weight();
  }
}


// Compute the source term for a single element
template <class LocalView, class LocalSourceTerm>
void getVolumeTerm( const LocalView& localView,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    LocalSourceTerm&& localSourceTerm)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();

  // Set all entries to zero
  localRhs.resize(localFiniteElement.size());
  localRhs = 0;

  // A quadrature rule
  int order = dim*localFiniteElement.localBasis().order();
  const auto& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad)
  {
    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPoint.position());

    double functionValue = localSourceTerm(quadPoint.position());

    // Evaluate all shape function values at this point
    std::vector<FieldVector<double,1> > shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPoint.position(), shapeFunctionValues);

    // Actually compute the vector entries
    for (size_t i=0; i<localRhs.size(); i++)
      localRhs[i] += shapeFunctionValues[i] * functionValue * quadPoint.weight() * integrationElement;
  }

}

// Get the occupation pattern of the stiffness matrix
template <class FEBasis>
void getOccupationPattern(const FEBasis& feBasis, MatrixIndexSet& nb)
{
  // Total number of grid vertices
  auto n = feBasis.size();

  nb.resize(n, n);

  // A view on the FE basis on a single element
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.localIndexSet();

  // Loop over all leaf elements
  for(const auto& e : elements(feBasis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(e);
    localIndexSet.bind(localView);

    for (size_t i=0; i<localView.tree().size(); i++) {

      for (size_t j=0; j<localView.tree().size(); j++) {

        auto iIdx = localIndexSet.index(i)[0];
        auto jIdx = localIndexSet.index(j)[0];

        // Add a nonzero entry to the matrix
        nb.add(iIdx, jIdx);

      }

    }

    // Now let's get the off-diagonal element stiffness matrix
    for (auto&& is : intersections(localView.globalBasis().gridView(), e))
    {
      if (!is.neighbor())
        continue;

      // Get a local view and local index set for the element on the other side of the intersection
      auto outsideLocalView = feBasis.localView();
      auto outsideLocalIndexSet = feBasis.localIndexSet();
      outsideLocalView.bind(is.outside());
      outsideLocalIndexSet.bind(outsideLocalView);

      // Add element stiffness matrix onto the global stiffness matrix
      for (size_t i=0; i<localView.tree().size(); i++)
      {
        // The global index of the i-th local degree of freedom of the element 'e'
        auto row = localIndexSet.index(i)[0];

        for (size_t j=0; j<outsideLocalView.tree().size(); j++ )
        {
          // The global index of the j-th local degree of freedom
          // of the element on the other side of the intersection
          auto col = outsideLocalIndexSet.index(j)[0];
          nb.add(row,col);
        }

      }

    }

  }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class FEBasis, class VelocityField, class ReactionCoefficient, class SourceTerm>
void assembleStiffnessMatrix(const FEBasis& feBasis,
                           BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                           BlockVector<FieldVector<double,1> >& rhs,
                           VelocityField&& velocityField,
                           ReactionCoefficient&& reactionCoefficient,
                           SourceTerm&& sourceTerm)
{
  // Get the grid view from the finite element basis
  typedef typename FEBasis::GridView GridView;
  GridView gridView = feBasis.gridView();

  // Make 'local' versions of the coefficient functions
  // 'local' means that you can bind them to grid elements, and evaluate them
  // in local coordinates of that element.
  auto localVelocityField       = localFunction(Functions::makeGridViewFunction(velocityField, gridView));
  auto localReactionCoefficient = localFunction(Functions::makeGridViewFunction(reactionCoefficient, gridView));
  auto localSourceTerm          = localFunction(Functions::makeGridViewFunction(sourceTerm, gridView));

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  getOccupationPattern(feBasis, occupationPattern);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(feBasis.size());

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // A view on the FE basis on a single element
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.localIndexSet();

  // A loop over all elements of the grid
  for(const auto& element : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    localVelocityField.bind(element);
    localReactionCoefficient.bind(element);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix, localVelocityField, localReactionCoefficient);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++) {

      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localIndexSet.index(i)[0];

      for (size_t j=0; j<elementMatrix.M(); j++ ) {

        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localIndexSet.index(j)[0];
        matrix[row][col] += elementMatrix[i][j];

      }

    }

    // Now let's get the off-diagonal element stiffness matrix
    for (auto&& is : intersections(gridView, element))
    {
      if (!is.neighbor())
        continue;

      // Get a local view and local index set for the element on the other side of the intersection
      auto outsideLocalView = feBasis.localView();
      auto outsideLocalIndexSet = feBasis.localIndexSet();
      outsideLocalView.bind(is.outside());
      outsideLocalIndexSet.bind(outsideLocalView);

      getOffDiagonalLocalMatrix(is, localView, outsideLocalView, elementMatrix, localVelocityField);

      // Add element stiffness matrix onto the global stiffness matrix
      for (size_t i=0; i<elementMatrix.N(); i++) {

        // The global index of the i-th local degree of freedom of the element 'e'
        auto row = localIndexSet.index(i)[0];

        for (size_t j=0; j<elementMatrix.M(); j++ ) {

          // The global index of the j-th local degree of freedom
          // of the element on the other side of the intersection
          auto col = outsideLocalIndexSet.index(j)[0];
          matrix[row][col] += elementMatrix[i][j];

        }
      }
    }

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    localSourceTerm.bind(element);
    getVolumeTerm(localView, localRhs, localSourceTerm);

    for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'e'
      auto row = localIndexSet.index(i)[0];
      rhs[row] += localRhs[i];

    }

  }

}


int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l = {1.0, 1.0};
  std::array<int,dim> elements = {10, 10};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  // Second-order Lagrange DG space
  typedef Functions::LagrangeDGBasis<GridView,2> FEBasis;
  FEBasis feBasis(gridView);

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  using Domain = GridType::Codim<0>::Geometry::GlobalCoordinate;  // FieldVector<float,dim>

  auto sourceTerm          = [] (const Domain& x) { return 10;};
  auto velocityField       = [] (const Domain& x) -> Domain { return {1,1};};
  auto reactionCoefficient = [] (const Domain& x) { return 10;};

  assembleStiffnessMatrix(feBasis, stiffnessMatrix, rhs, velocityField, reactionCoefficient, sourceTerm);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasis.size());
  x = 0;

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  SeqILU0<MatrixType,VectorType,VectorType> preconditioner(stiffnessMatrix,1.0);

//  Richardson<VectorType,VectorType> preconditioner(1.0);

  // Preconditioned conjugate-gradient solver
  RestartedGMResSolver<VectorType> solver(op,
                                          preconditioner,
                                          1e-10,  // desired residual reduction factor
                                          500,     // number of iterations between restarts
                                          500,   // maximum number of iterations
                                          2);    // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  solver.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  auto xFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(feBasis, Dune::TypeTree::hybridTreePath(), x);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(xFunction, VTK::FieldInfo("x", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("advection-reaction-dg");

}
// Error handling
catch (Exception e)
{
  std::cout << e << std::endl;
}
