// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>
#include <cmath>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/bdmbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#define RT0
//#define BDM1
//#define DIM2

using namespace Dune;

// Compute the stiffness matrix for a single element
template <class LocalView, class MatrixType>
void getLocalMatrix(const LocalView& localView,
                    MatrixType& elementMatrix)
{
  // Get the grid element from the local FE basis view (use test space)
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim   = Element::dimension;
  auto geometry = element.geometry();

  // Set all matrix entries to zero
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // Get set of shape functions for this element
  using namespace Dune::TypeTree::Indices;
  const auto& fluxLocalFiniteElement     = localView.tree().child(_0).finiteElement();
  const auto& pressureLocalFiniteElement = localView.tree().child(_1).finiteElement();

  // Get a quadrature rule
  int order = std::max(2*(dim*fluxLocalFiniteElement.localBasis().order()-1), 2*(dim*pressureLocalFiniteElement.localBasis().order()));
  const auto& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad) {

    // Position of the current quadrature point in the reference element
    const auto quadPos = quadPoint.position();

    // Values of the pressure shape functions
    std::vector<FieldVector<double,1> > pressureValues;
    pressureLocalFiniteElement.localBasis().evaluateFunction(quadPos, pressureValues);

    // Values of the flux shape functions
    std::vector<FieldVector<double,dim> > fluxValues;
    fluxLocalFiniteElement.localBasis().evaluateFunction(quadPos, fluxValues);

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    // Gradients of the flux shape functions on the reference element
    std::vector<FieldMatrix<double,dim,dim> > fluxReferenceGradients;
    fluxLocalFiniteElement.localBasis().evaluateJacobian(quadPos, fluxReferenceGradients);

    // Compute the flux shape function gradients on the real element
    std::vector<FieldMatrix<double,dim,dim> > fluxGradients(fluxReferenceGradients.size());
    for (size_t i=0; i<fluxGradients.size(); i++)
      for (size_t j=0; j<dim; j++)
        jacobian.mv(fluxReferenceGradients[i][j], fluxGradients[i][j]);

    /////////////////////////////
    // Flux--flux coupling
    /////////////////////////////

    for (size_t i=0; i<fluxLocalFiniteElement.size(); i++) {
      size_t row = localView.tree().child(_0).localIndex(i);
      for (size_t j=0; j<fluxLocalFiniteElement.size(); j++) {
          size_t col = localView.tree().child(_0).localIndex(j);
          elementMatrix[row][col] += (fluxValues[i] * fluxValues[j]) * quadPoint.weight() * integrationElement;
      }
    }

    /////////////////////////////
    // Flux--pressure coupling
    /////////////////////////////

    for (size_t i=0; i<fluxLocalFiniteElement.size(); i++) {
      size_t fluxIndex     = localView.tree().child(_0).localIndex(i);

      // Pre-compute divergence
      double fluxDivergence = 0.;
      for (size_t k=0; k<dim; k++)
        fluxDivergence += fluxGradients[i][k][k];

      for (size_t j=0; j<pressureLocalFiniteElement.size(); j++) {
        size_t pressureIndex = localView.tree().child(_1).localIndex(j);

        // Pre-compute matrix contribution
        double tmp = - (pressureValues[j] * fluxDivergence) * quadPoint.weight() * integrationElement;

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
  using namespace Dune::TypeTree::Indices;
  const auto& fluxLocalFiniteElement     = localView.tree().child(_0).finiteElement();
  const auto& pressureLocalFiniteElement = localView.tree().child(_1).finiteElement();

  // A quadrature rule
  int order = std::max(dim*fluxLocalFiniteElement.localBasis().order(), dim*pressureLocalFiniteElement.localBasis().order());
  const auto& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad) {

    // Position of the current quadrature point in the reference element
    const auto& quadPos = quadPoint.position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPos);

    double functionValue = localVolumeTerm(quadPos);

    // Evaluate all shape function values at this point
    std::vector<FieldVector<double,1> > pressureValues;
    pressureLocalFiniteElement.localBasis().evaluateFunction(quadPos, pressureValues);

    // Actually compute the vector entries
    for (size_t j=0; j<pressureLocalFiniteElement.size(); j++) {
      size_t pressureIndex = localView.tree().child(_1).localIndex(j);
      localRhs[pressureIndex] += (-1) * pressureValues[j] * functionValue * quadPoint.weight() * integrationElement;

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
  auto localIndexSet = basis.localIndexSet();

  // Loop over all leaf elements
  for(const auto& element : elements(basis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    // Add element stiffness matrix onto global stiffness matrix
    for (size_t i=0; i<localIndexSet.size(); i++) {

      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<localIndexSet.size(); ++j) {

        // The global index set of the j-th local degree of freedom of the element 'e'
        auto col = localIndexSet.index(j);

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
  typedef typename Basis::GridView GridView;
  GridView gridView = basis.gridView();

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
  auto localIndexSet = basis.localIndexSet();

  // A loop over all elements of the grid
  for(const auto& element : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++) {

      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<elementMatrix.M(); j++ ) {

        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localIndexSet.index(j);
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
  typedef typename Basis::GridView GridView;
  GridView gridView = basis.gridView();

  auto localVolumeTerm = localFunction(Functions::makeGridViewFunction(volumeTerm, gridView));

  // set rhs to correct length -- the total number of basis vectors in the basis
  typedef Dune::Functions::HierarchicVectorWrapper<VectorType, double> HierarchicVectorView;
  HierarchicVectorView(rhs).resize(basis);

  // Set all entries to zero
  rhs = 0;

  // A view on the FE basis on a single element
  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();

  // A loop over all elements of the grid
  for(const auto& element : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    localVolumeTerm.bind(element);
    getVolumeTerm(localView, localRhs, localVolumeTerm);

    for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'e'
      auto row = localIndexSet.index(i);
      rhs[row[0]][row[1]] += localRhs[i];

    }

  }

}

int main (int argc, char *argv[]) try
{
#ifndef RT0
#ifndef BDM1
  DUNE_THROW(Dune::NotImplemented, "Choose RT0 or BDM1.");
#endif
#endif

  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

#ifdef DIM2
  const int dim = 2;
  std::array<int,dim> elements = {20, 20};
#else
  const int dim = 3;
  std::array<int,dim> elements = {20, 20, 20};
#endif
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisTags;
  using namespace Functions::BasisBuilder;

  auto basis = makeBasis(
    gridView,
    composite(
#ifdef RT0
      rt<0, GeometryType::BasicType::cube>(),
      pq<0>()
#endif
#ifdef BDM1
      bdm<1, GeometryType::BasicType::cube>(),
      pq<0>()
#endif
    ));


  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  typedef Matrix<BCRSMatrix<FieldMatrix<double,1,1> > > MatrixType;
  typedef BlockVector<BlockVector<FieldVector<double,1> > > VectorType;
  typedef Dune::Functions::HierarchicVectorWrapper<VectorType, double> HierarchicVectorView;

  MatrixType stiffnessMatrix;
  VectorType rhs;

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;

  assembleMixedPoissonMatrix(basis, stiffnessMatrix);

  auto rightHandSide = [] (const Domain& x) { return 2.;};
  assembleMixedPoissonRhs(basis, rhs, rightHandSide);

  auto topFluxBC = [] (const Domain& x) { return 0.0; };
  auto lowerFluxBC = [] (const Domain& x) { return 0.0; };

  using namespace TypeTree::Indices;
  using BitVectorType = BlockVector<BlockVector<FieldVector<char,1> > >;
  using HierarchicBitVectorView = Functions::HierarchicVectorWrapper<BitVectorType, char>;

  BitVectorType isTopBoundary;
  BitVectorType isLowerBoundary;

#ifdef RT0

  // Use the same way as for PQk functions to mark boundaries.

  // Mark top boundary.
  auto topBoundaryIndicator = [&l] (Domain x) {
        bool isBoundary = x[dim-1] > l[dim-1] - 1e-8;
        return isBoundary;
  };

  // Mark top boundary.
  auto lowerBoundaryIndicator = [&l] (Domain x) {
        bool isBoundary = x[dim-1] < 1e-8;
        return isBoundary;
  };

  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicBitVectorView(isTopBoundary), topBoundaryIndicator);
  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicBitVectorView(isLowerBoundary), lowerBoundaryIndicator);

#else

  // Use a messy way of defining a boundary. Using the RT0-way, only the RT0 basis DOF are triggered,
  // as they suffice to interpolate a constant function on an edge. In particular BDM1 DOF are not triggered,
  // which are linear (but non-constant) on edges.

  // Mark top boundary.
  auto topBoundaryIndicator = [&l] (Domain x) {
       double isBoundary = x[dim-1] > l[dim-1] - 1e-8 ? x[0] : 0.0;
        return isBoundary;
  };

  // Mark top boundary.
  auto lowerBoundaryIndicator = [&l] (Domain x) {
        double isBoundary = x[dim-1] < 1e-8 ? x[0] : 0.0;
        return isBoundary;
  };

  VectorType isTopBoundaryTmp, isLowerBoundaryTmp;

  // Use double-valued interpolation and transfer to char-valued vectors.
  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(isTopBoundaryTmp), topBoundaryIndicator);
  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(isLowerBoundaryTmp), lowerBoundaryIndicator);
  HierarchicBitVectorView(isTopBoundary).resize(basis);
  HierarchicBitVectorView(isLowerBoundary).resize(basis);
  isTopBoundary = 0;
  isLowerBoundary = 0;
  for (size_t i=0; i<isTopBoundaryTmp[0].size(); i++){
    isTopBoundary[0][i] = isTopBoundaryTmp[0][i]!=0 ? 1: 0;
    isLowerBoundary[0][i] = isLowerBoundaryTmp[0][i]!=0 ? 1: 0;
  }
#endif

  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(rhs), topFluxBC, HierarchicBitVectorView(isTopBoundary));
  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(rhs), lowerFluxBC, HierarchicBitVectorView(isLowerBoundary));

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<stiffnessMatrix[0][0].N(); i++) {

    if (isTopBoundary[0][i] or isLowerBoundary[0][i]) {

      // Lower right matrix block
      auto cIt    = stiffnessMatrix[0][0][i].begin();
      auto cEndIt = stiffnessMatrix[0][0][i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt){
        *cIt = (i==cIt.index()) ? 1. : 0.;
      }

      // Lower left matrix block
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

  // Preconditioned conjugate-gradient solver
  RestartedGMResSolver<VectorType> solver(op,
                                          preconditioner,
                                          1e-6,  // desired residual reduction factor
                                          100,     // number of iterations between restarts
                                          1000,   // maximum number of iterations
                                          2);    // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  solver.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  using FluxRange = FieldVector<double,dim>;
  using PressureRange = double;

  auto fluxFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<FluxRange>(basis, TypeTree::hybridTreePath(_0), HierarchicVectorView(x));

  auto pressureFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<PressureRange>(basis, TypeTree::hybridTreePath(_1), HierarchicVectorView(x));

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(fluxFunction, VTK::FieldInfo("flux", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.addCellData(pressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("poisson-mfem-result");
 }
// Error handling
 catch (Exception e) {
    std::cout << e << std::endl;
 }
