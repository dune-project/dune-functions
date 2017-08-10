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
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#define DIM2 // Use a two-dimensional test, otherwise three-dimensional

using namespace Dune;

// Compute the stiffness matrix for a single element
template <class LocalView, class MatrixType>
void getLocalMatrix(const LocalView& localView,
                    MatrixType& elementMatrix)
{
  // Get the grid element from the local FE basis view (use test space)
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;
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
  for (const auto& quadPoint : quad)
  {
    // Position of the current quadrature point in the reference element
    const auto quadPos = quadPoint.position();

    // The transposed (inverse) Jacobian of the map from the reference element to the element
    const auto& jacobianTransposed = geometry.jacobianTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula and its inverse
    const auto integrationElement = geometry.integrationElement(quadPos);
    const auto invIntegrationElement = 1.0 / integrationElement;

    ///////////////////////////////////////////////////////////////////////////
    // Shape functions - flux - [assume H(div) elements]
    ///////////////////////////////////////////////////////////////////////////

    // Values of the flux shape functions on the reference element
    std::vector<FieldVector<double,dim> > fluxReferenceValues(fluxLocalFiniteElement.size());
    fluxLocalFiniteElement.localBasis().evaluateFunction(quadPos, fluxReferenceValues);

    // Values of the shape functions on the current element
    std::vector<FieldVector<double,dim> > fluxValues(fluxLocalFiniteElement.size());
    for (size_t i=0; i<fluxValues.size(); i++)
    {
      jacobianTransposed.mtv(fluxReferenceValues[i], fluxValues[i]);
      fluxValues[i] *= invIntegrationElement;
    }

    // Gradients of the flux shape functions on the reference element
    std::vector<FieldMatrix<double,dim,dim> > fluxReferenceGradients(fluxLocalFiniteElement.size());
    fluxLocalFiniteElement.localBasis().evaluateJacobian(quadPos, fluxReferenceGradients);

    // Compute the flux shape function divergence on the real element
    std::vector<double> fluxDivergence(fluxValues.size(), 0.0);
    for (size_t i=0; i<fluxValues.size(); i++)
    {
      for (size_t j=0; j<dim; j++)
        fluxDivergence[i] += fluxReferenceGradients[i][j][j];
      fluxDivergence[i] *= invIntegrationElement;
    }

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
        double tmp = - (fluxDivergence[i] * pressureValues[j]) * quadPoint.weight() * integrationElement;

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
  auto localIndexSet = basis.localIndexSet();

  // Loop over all leaf elements
  for(const auto& element : elements(basis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);

    // Add element stiffness matrix onto global stiffness matrix
    for (size_t i=0; i<localIndexSet.size(); i++)
    {
      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<localIndexSet.size(); ++j)
      {
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
    for (size_t i=0; i<elementMatrix.N(); i++)
    {
      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<elementMatrix.M(); j++ )
      {
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

    for (size_t i=0; i<localRhs.size(); i++)
    {
      // The global index of the i-th vertex of the element 'e'
      auto row = localIndexSet.index(i);
      rhs[row[0]][row[1]] += localRhs[i];
    }
  }
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
  std::array<int,dim> elements = {50, 50};
#else
  const int dim = 3;
  std::array<int,dim> elements = {10, 10, 10};
#endif
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisBuilder;

  const int k = 0;

  auto basis = makeBasis(
    gridView,
    composite(
      rt<k, GeometryType::BasicType::cube>(),
      pq<k>()
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

  auto rightHandSide = [] (const Domain& x) { return 10; };
  assembleMixedPoissonRhs(basis, rhs, rightHandSide);

  auto pi = std::acos(-1.0);
  auto topFluxBC   = [&pi] (const Domain& x) { return -0.05 * (1. - x[0]) * std::sin(2.*pi*x[0]); };
  auto lowerFluxBC = [&pi] (const Domain& x) { return  0.05 * (1. - x[0]) * std::sin(2.*pi*x[0]); };

  using namespace TypeTree::Indices;
  using BitVectorType = BlockVector<BlockVector<FieldVector<char,1> > >;
  using HierarchicBitVectorView = Functions::HierarchicVectorWrapper<BitVectorType, char>;

  BitVectorType isTopBoundary;
  BitVectorType isLowerBoundary;

  // Use a messy way of defining a boundary. Using the same way as for PQk elements, only the RT0 basis DOF are triggered,
  // as they suffice to interpolate a constant function on an edge. In particular RT1, BDM1 DOF are not triggered,
  // which are linear (but non-constant) on edges.

  // Mark top boundary.
  auto topBoundaryIndicator = [&l] (Domain x)
  {
#ifdef DIM2
    double isBoundary = x[dim-1] > l[dim-1] - 1e-8 ? x[0]: 0.0;
#else
    double isBoundary = x[dim-1] > l[dim-1] - 1e-8 ? x[0] * x[1]: 0.0;
#endif
    return isBoundary;
  };

  // Mark lower boundary.
  auto lowerBoundaryIndicator = [&l] (Domain x)
  {
#ifdef DIM2
    double isBoundary = x[dim-1] < 1e-8 ? x[0]: 0.0;
#else
    double isBoundary = x[dim-1] < 1e-8 ? x[0] * x[1]: 0.0;
#endif
    return isBoundary;
  };

  VectorType isTopBoundaryTmp, isLowerBoundaryTmp;

  // Use double-valued interpolation and transfer to char-valued vectors.
  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(isTopBoundaryTmp),   topBoundaryIndicator);
  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(isLowerBoundaryTmp), lowerBoundaryIndicator);
  HierarchicBitVectorView(isTopBoundary).resize(basis);
  HierarchicBitVectorView(isLowerBoundary).resize(basis);
  isTopBoundary = 0;
  isLowerBoundary = 0;
  for (size_t i=0; i<isTopBoundaryTmp[0].size(); i++)
  {
    isTopBoundary[  0][i] = isTopBoundaryTmp[  0][i]!=0 ? 1: 0;
    isLowerBoundary[0][i] = isLowerBoundaryTmp[0][i]!=0 ? 1: 0;
  }

  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(rhs), topFluxBC,   HierarchicBitVectorView(isTopBoundary));
  interpolate(basis, Dune::TypeTree::hybridTreePath(_0), HierarchicVectorView(rhs), lowerFluxBC, HierarchicBitVectorView(isLowerBoundary));

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<stiffnessMatrix[0][0].N(); i++)
  {
    if (isTopBoundary[0][i] or isLowerBoundary[0][i])
    {
      // Lower right matrix block
      auto cIt    = stiffnessMatrix[0][0][i].begin();
      auto cEndIt = stiffnessMatrix[0][0][i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? 1. : 0.;

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

  auto fluxFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<FluxRange>(basis, TypeTree::hybridTreePath(_0), HierarchicVectorView(x));
  auto pressureFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<PressureRange>(basis, TypeTree::hybridTreePath(_1), HierarchicVectorView(x));

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,0);
  vtkWriter.addVertexData(fluxFunction, VTK::FieldInfo("flux", VTK::FieldInfo::Type::vector, dim));
  if (k==0)
    vtkWriter.addCellData(pressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  else
    vtkWriter.addVertexData(pressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("poisson-mfem-result");
}
