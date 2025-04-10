// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <vector>
#include <cmath>

#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/facenormalgridfunction.hh>
#include <dune/functions/gridfunctions/composedgridfunction.hh>

#define DIM2 // Use a two-dimensional test, otherwise three-dimensional

using namespace Dune;

// Compute the stiffness matrix for a single element
template <class LocalView, class MatrixType>
void getLocalMatrix(const LocalView& localView,
                    MatrixType& elementMatrix)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Set all matrix entries to zero
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // Get set of shape functions for this element
  using namespace Dune::Indices;
  const auto& fluxLocalFiniteElement     = localView.tree().child(_0).finiteElement();
  const auto& pressureLocalFiniteElement = localView.tree().child(_1).finiteElement();

  // Get a quadrature rule
  int fluxOrder = dim*fluxLocalFiniteElement.localBasis().order();
  int pressureOrder = dim*pressureLocalFiniteElement.localBasis().order();
  int order = std::max(2*fluxOrder, (fluxOrder-1)+pressureOrder);
  const auto& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad)
  {
    // Position of the current quadrature point in the reference element
    const auto quadPos = quadPoint.position();

    // The inverse Jacobian of the map from the reference element to the element
    const auto jacobianInverse = geometry.jacobianInverse(quadPos);

    // The multiplicative factor in the integral transformation formula
    const auto integrationElement = geometry.integrationElement(quadPos);

    ///////////////////////////////////////////////////////////////////////////
    // Shape functions - flux
    ///////////////////////////////////////////////////////////////////////////

    // Values of the flux shape functions on the current element
    std::vector<FieldVector<double,dim> > fluxValues(fluxLocalFiniteElement.size());
    fluxLocalFiniteElement.localBasis().evaluateFunction(quadPos, fluxValues);

    // Gradients of the flux shape function gradients on the reference element
    std::vector<FieldMatrix<double,dim,dim> > fluxReferenceJacobians(fluxLocalFiniteElement.size());
    fluxLocalFiniteElement.localBasis().evaluateJacobian(quadPos, fluxReferenceJacobians);

    // Helper function to compute the trace of a matrix
    auto trace = [](const auto& matrix) {
      double r=0;
      for (size_t j=0; j<matrix.N(); j++)
        r += matrix[j][j];
      return r;
    };

    // Domain transformation of Jacobians and computation of div = trace(Jacobian)
    std::vector<double> fluxDivergence(fluxValues.size(), 0.0);
    for (size_t i=0; i<fluxReferenceJacobians.size(); i++)
      fluxDivergence[i] = trace(fluxReferenceJacobians[i] * jacobianInverse);

    ///////////////////////////////////////////////////////////////////////////
    // Shape functions - pressure
    ///////////////////////////////////////////////////////////////////////////

    // Values of the pressure shape functions
    std::vector<FieldVector<double,1> > pressureValues(pressureLocalFiniteElement.size());
    pressureLocalFiniteElement.localBasis().evaluateFunction(quadPos, pressureValues);

    /////////////////////////////
    // Flux--flux coupling
    /////////////////////////////

    for (size_t i=0; i<fluxLocalFiniteElement.size(); i++)
    {
      size_t row = localView.tree().child(_0).localIndex(i);
      for (size_t j=0; j<fluxLocalFiniteElement.size(); j++)
      {
        size_t col = localView.tree().child(_0).localIndex(j);
        elementMatrix[row][col] += (fluxValues[i] * fluxValues[j]) * quadPoint.weight() * integrationElement;
      }
    }

    /////////////////////////////
    // Flux--pressure coupling
    /////////////////////////////

    for (size_t i=0; i<fluxLocalFiniteElement.size(); i++)
    {
      size_t fluxIndex     = localView.tree().child(_0).localIndex(i);
      for (size_t j=0; j<pressureLocalFiniteElement.size(); j++)
      {
        size_t pressureIndex = localView.tree().child(_1).localIndex(j);

        // Pre-compute matrix contribution
        double tmp = (fluxDivergence[i] * pressureValues[j]) * quadPoint.weight() * integrationElement;

        elementMatrix[fluxIndex][pressureIndex] += tmp;
        elementMatrix[pressureIndex][fluxIndex] += tmp;
      }
    }
  }
}


// Compute the source term for a single element
template <class LocalView, class LocalVolumeTerm>
void getVolumeTerm( const LocalView& localView,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    LocalVolumeTerm&& localVolumeTerm)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;

  // Set all entries to zero
  localRhs.resize(localView.size());
  localRhs = 0;

  // Get set of shape functions for this element
  using namespace Dune::Indices;
  const auto& fluxLocalFiniteElement     = localView.tree().child(_0).finiteElement();
  const auto& pressureLocalFiniteElement = localView.tree().child(_1).finiteElement();

  // A quadrature rule
  int fluxOrder = dim*fluxLocalFiniteElement.localBasis().order();
  int pressureOrder = dim*pressureLocalFiniteElement.localBasis().order();
  int order = std::max(2*fluxOrder, 2*pressureOrder);
  const auto& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad)
  {
    // Position of the current quadrature point in the reference element
    const auto& quadPos = quadPoint.position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPos);

    double functionValue = localVolumeTerm(quadPos);

    // Evaluate all shape function values at this point
    std::vector<FieldVector<double,1> > pressureValues;
    pressureLocalFiniteElement.localBasis().evaluateFunction(quadPos, pressureValues);

    // Actually compute the vector entries
    for (size_t j=0; j<pressureLocalFiniteElement.size(); j++)
    {
      size_t pressureIndex = localView.tree().child(_1).localIndex(j);
      localRhs[pressureIndex] += - pressureValues[j] * functionValue * quadPoint.weight() * integrationElement;
    }
  }
}

// Get the occupation pattern of the stiffness matrix
template <class Basis>
void getOccupationPattern(const Basis& basis,
                          std::array<std::array<MatrixIndexSet,2>,2>& nb)
{
  // Set sizes of the 2x2 submatrices
  for (size_t i=0; i<2; ++i)
    for (size_t j=0; j<2; j++)
      nb[i][j].resize(basis.size({i}), basis.size({j}));

  // A view on the FE basis on a single element
  auto localView = basis.localView();

  // Loop over all leaf elements
  for(const auto& element : elements(basis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);

    // Add element stiffness matrix onto global stiffness matrix
    for (size_t i=0; i<localView.size(); i++)
    {
      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localView.index(i);

      for (size_t j=0; j<localView.size(); ++j)
      {
        // The global index set of the j-th local degree of freedom of the element 'e'
        auto col = localView.index(j);

        nb[row[0]][col[0]].add(row[1], col[1]);
      }
    }
  }
}

/** \brief Assemble the divergence stiffness matrix on the given grid view */
template <class Basis, class MatrixType>
void assembleMixedPoissonMatrix(const Basis& basis,
                                MatrixType& matrix)
{
  // Get the grid view from the finite element basis (use of test space arbitrary)
  auto gridView = basis.gridView();

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  std::array<std::array<MatrixIndexSet, 2>, 2> occupationPattern;
  getOccupationPattern(basis, occupationPattern);
  // ... and give it the occupation pattern we want.
  matrix.setSize(2,2);
  for (size_t i=0; i<2; i++)
    for (size_t j=0; j<2; j++)
      occupationPattern[i][j].exportIdx(matrix[i][j]);

  // Set all entries to zero
  matrix = 0;

  // A view on the FE basis on a single element
  auto localView = basis.localView();

  // A loop over all elements of the grid
  for(const auto& element : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++)
    {
      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localView.index(i);

      for (size_t j=0; j<elementMatrix.M(); j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localView.index(j);
        matrix[row[0]][col[0]][row[1]][col[1]] += elementMatrix[i][j];
      }
    }
  }
}

/** \brief Assemble the divergence stiffness matrix on the given grid view */
template <class Basis, class VectorType, class VolumeTerm>
void assembleMixedPoissonRhs(const Basis& basis,
                             VectorType& rhs,
                             VolumeTerm&& volumeTerm)
{
  // Get the grid view from the finite element basis (use of test space arbitrary)
  auto gridView = basis.gridView();

  auto localVolumeTerm = localFunction(Functions::makeGridViewFunction(volumeTerm, gridView));

  // set rhs to correct length -- the total number of basis vectors in the basis
  using Functions::istlVectorBackend;
  istlVectorBackend(rhs).resize(basis);

  // Set all entries to zero
  rhs = 0;

  // A view on the FE basis on a single element
  auto localView = basis.localView();

  // A loop over all elements of the grid
  for(const auto& element : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    localVolumeTerm.bind(element);
    getVolumeTerm(localView, localRhs, localVolumeTerm);

    for (size_t i=0; i<localRhs.size(); i++)
    {
      // The global index of the i-th vertex of the element 'e'
      auto row = localView.index(i);
      rhs[row[0]][row[1]] += localRhs[i];
    }
  }
}



// Mark all DOFs associated to entities for which
// the boundary intersections center is marked
// by the given indicator functions.
template<class Basis, class Vector, class Indicator>
void markBoundaryDOFsByIndicator(const Basis& basis, Vector& vector, const Indicator& indicator)
{
  auto vectorBackend = Dune::Functions::istlVectorBackend(vector);
  Dune::Functions::forEachBoundaryDOF(basis, [&] (auto&& localIndex, const auto& localView, const auto& intersection) {
    if (indicator(intersection.geometry().center())>1e-8)
      vectorBackend[localView.index(localIndex)] = true;
  });
}



int main (int argc, char *argv[])
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

#ifdef DIM2
  const int dim = 2;
  std::array<int,dim> elements = {{50, 50}};
#else
  const int dim = 3;
  std::array<int,dim> elements = {{10, 10, 10}};
#endif
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  GridType grid(l,elements);

  auto gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisFactory;
  using namespace Indices;

  const int k = 0;

  auto basis = makeBasis(
    gridView,
    composite(
      raviartThomas<k>(),
      lagrange<k>(),
      blockedLexicographic()
    ));

  auto fluxBasis = Functions::subspaceBasis(basis, _0);
  auto pressureBasis = Functions::subspaceBasis(basis, _1);


  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  using MatrixType = Matrix<BCRSMatrix<double> >;
  using VectorType = BlockVector<BlockVector<double> >;
  using BitVectorType = BlockVector<BlockVector<char> >;

  using Functions::istlVectorBackend;

  MatrixType stiffnessMatrix;
  VectorType rhs;

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;

  assembleMixedPoissonMatrix(basis, stiffnessMatrix);

  auto rightHandSide = [] (const Domain& x) { return 2; };

  assembleMixedPoissonRhs(basis, rhs, rightHandSide);

  // This marks the top and bottom boundary of the domain
  auto fluxDirichletIndicator = [&l] (const auto& x) {
    return ((x[dim-1] > l[dim-1] - 1e-8) or (x[dim-1] < 1e-8));
  };

  auto coordinate = Dune::Functions::makeAnalyticGridViewFunction([](const auto& x) { return x; }, gridView);
  auto normal = Dune::Functions::FaceNormalGridFunction(gridView);
  auto fluxDirichletValues = Dune::Functions::makeComposedGridFunction(
    [pi = std::acos(-1.0)](const auto& x, const auto& normal) {
      return std::sin(2.*pi*x[0]) * normal;
    },
    coordinate,
    normal);

  // Mark all DOFs located in a boundary intersection marked
  // by the fluxDirichletIndicator function. If the flux
  // ansatz space also contains tangential components, this
  // approach will fail, because those are also marked.
  // For Raviart-Thomas this does not happen.
  auto isDirichlet = BitVectorType();
  istlVectorBackend(isDirichlet).resize(basis);
  isDirichlet = false;
  markBoundaryDOFsByIndicator(fluxBasis, isDirichlet, fluxDirichletIndicator);

  // Interpolate flux Dirichlet values
  interpolate(fluxBasis, rhs, fluxDirichletValues, isDirichlet);

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<stiffnessMatrix[0][0].N(); i++)
  {
    if (isDirichlet[0][i])
    {
      // Upper left matrix block
      for (auto&& [entry, idx] : sparseRange(stiffnessMatrix[0][0][i]))
        entry = (i==idx) ? 1.0 : 0.0;

      // Upper right matrix block
      for (auto&& entry: stiffnessMatrix[0][1][i])
        entry = 0.0;
    }
  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Start from the rhs vector; that way the Dirichlet entries are already correct
  VectorType x = rhs;

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Fancy (but only) way to not have a preconditioner at all
  Richardson<VectorType,VectorType> preconditioner(1.0);

  // Preconditioned GMRES / BiCGSTAB solver
  //RestartedGMResSolver<VectorType> solver (op, preconditioner, 1e-6, 1000, 10000, 2);
  BiCGSTABSolver<VectorType> solver(op, preconditioner, 1e-6, 4000, 2);

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  solver.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  using FluxRange = FieldVector<double,dim>;
  using PressureRange = double;

  auto fluxFunction = Functions::makeDiscreteGlobalBasisFunction<FluxRange>(fluxBasis, x);
  auto pressureFunction = Functions::makeDiscreteGlobalBasisFunction<PressureRange>(pressureBasis, x);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //////////////////////////////////////////////////////////////////////////////////////////////
  auto vtkWriter = SubsamplingVTKWriter(gridView, Dune::refinementLevels(0));
  vtkWriter.addVertexData(fluxFunction, VTK::FieldInfo("flux", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.addVertexData(pressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("poisson-mfem-result");
}
