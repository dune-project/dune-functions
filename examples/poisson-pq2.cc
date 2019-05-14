// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>
#include <cmath>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

using namespace Dune;

// Compute the stiffness matrix for a single element
template <class LocalView, class MatrixType>
void getLocalMatrix( const LocalView& localView, MatrixType& elementMatrix)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();

  // Set all matrix entries to zero
  elementMatrix.setSize(localFiniteElement.localBasis().size(),localFiniteElement.localBasis().size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // Get a quadrature rule
  int order = 2*(dim*localFiniteElement.localBasis().order()-1);
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (size_t pt=0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos, referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
    for (size_t i=0; i<gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);

    // Compute the actual matrix entries
    for (size_t i=0; i<elementMatrix.N(); i++)
      for (size_t j=0; j<elementMatrix.M(); j++ )
        elementMatrix[i][j] += ( gradients[i] * gradients[j] ) * quad[pt].weight() * integrationElement;

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

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();

  // Set all entries to zero
  localRhs.resize(localFiniteElement.localBasis().size());
  localRhs = 0;

  // A quadrature rule
  int order = dim*localFiniteElement.localBasis().order();
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for ( size_t pt=0; pt < quad.size(); pt++ ) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPos);

    double functionValue = localVolumeTerm(quadPos);

    // Evaluate all shape function values at this point
    std::vector<FieldVector<double,1> > shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    // Actually compute the vector entries
    for (size_t i=0; i<localRhs.size(); i++)
      localRhs[i] += shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;

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

  // Loop over all leaf elements
  for(const auto& e : elements(feBasis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(e);

    // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
    for (size_t i=0; i<localView.tree().size(); i++) {

      for (size_t j=0; j<localView.tree().size(); j++) {

        auto iIdx = localView.index(i);
        auto jIdx = localView.index(j);

        // Add a nonzero entry to the matrix
        nb.add(iIdx, jIdx);

      }

    }

  }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class FEBasis, class VolumeTerm>
void assembleLaplaceMatrix(const FEBasis& feBasis,
                           BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                           BlockVector<FieldVector<double,1> >& rhs,
                           VolumeTerm&& volumeTerm)
{
  // Get the grid view from the finite element basis
  typedef typename FEBasis::GridView GridView;
  GridView gridView = feBasis.gridView();

  auto localVolumeTerm = localFunction(Functions::makeGridViewFunction(volumeTerm, gridView));




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

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView))
  {

    // Bind the local FE basis view to the current element
    localView.bind(e);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++) {

      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localView.index(i);

      for (size_t j=0; j<elementMatrix.M(); j++ ) {

        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localView.index(j);
        matrix[row][col] += elementMatrix[i][j];

      }

    }

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    localVolumeTerm.bind(e);
    getVolumeTerm(localView, localRhs, localVolumeTerm);

    for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'e'
      auto row = localView.index(i);
      rhs[row] += localRhs[i];

    }

  }

}


// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
// Since interpolating into a vector<bool> is currently not supported,
// we use a vector<char> which, in contrast to vector<bool>
// is a real container.
template <class FEBasis>
void boundaryTreatment (const FEBasis& feBasis, std::vector<char>& dirichletNodes )
{
  dirichletNodes.clear();
  dirichletNodes.resize(feBasis.size(), false);

  Functions::forEachBoundaryDOF(feBasis, [&] (auto&& index) {
    dirichletNodes[index] = true;
  });
}

template<int dim>
auto createUniformCubeGrid()
{
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> l(1.0);
  std::array<int,dim> elements = {{2, 2}};
  return std::make_unique<Grid>(l, elements);
}

#if HAVE_UG
auto createMixedGrid()
{
  using Grid = Dune::UGGrid<2>;
  Dune::GridFactory<Grid> factory;
  for(unsigned int k : Dune::range(9))
    factory.insertVertex({0.5*(k%3), 0.5*(k/3)});
  factory.insertElement(Dune::GeometryTypes::cube(2), {0, 1, 3, 4});
  factory.insertElement(Dune::GeometryTypes::cube(2), {1, 2, 4, 5});
  factory.insertElement(Dune::GeometryTypes::simplex(2), {3, 4, 6});
  factory.insertElement(Dune::GeometryTypes::simplex(2), {4, 7, 6});
  factory.insertElement(Dune::GeometryTypes::simplex(2), {4, 5, 7});
  factory.insertElement(Dune::GeometryTypes::simplex(2), {5, 8, 7});
  return std::unique_ptr<Grid>{factory.createGrid()};
}
#endif

int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

#if HAVE_UG
  auto gridPtr = createMixedGrid();
#else
  auto gridPtr = createUniformCubeGrid<2>();
#endif

  auto& grid = *gridPtr;
  grid.globalRefine(4);

  auto gridView = grid.leafGridView();
  using GridView = decltype(gridView);

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::LagrangeBasis<GridView,2> FEBasis;
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

  std::cout << "Number of DOFs is " << feBasis.dimension() << std::endl;

  auto rightHandSide = [] (const auto& x) { return 10;};

  Dune::Timer timer;
  assembleLaplaceMatrix(feBasis, stiffnessMatrix, rhs, rightHandSide);
  std::cout << "Assembling the problem took " << timer.elapsed() << "s" << std::endl;

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasis.size());
  x = 0;

  // Determine Dirichlet dofs
  std::vector<char> dirichletNodes;
  boundaryTreatment(feBasis, dirichletNodes);



  // Don't trust on non-standard M_PI.
  auto pi = std::acos(-1.0);
  auto dirichletValueFunction = [pi](const auto& x){ return std::sin(2*pi*x[0]); };

  // Interpolate dirichlet values at the boundary nodes
  interpolate(feBasis, x, dirichletValueFunction, dirichletNodes);

  //////////////////////////////////////////////////////
  //   Incorporate dirichlet values in a symmetric way
  //////////////////////////////////////////////////////

  // Compute residual for non-homogeneous Dirichlet values
  // stored in x. Since x is zero for non-Dirichlet DOFs,
  // we can simply multiply by the matrix.
  stiffnessMatrix.mmv(x, rhs);

  // Change Dirichlet matrix rows and columns to the identity.
  for (size_t i=0; i<stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i]) {
      rhs[i] = x[i];
      auto cIt    = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? 1 : 0;
    }
    else {
      auto cIt    = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        if (dirichletNodes[cIt.index()])
          *cIt = 0;
    }
  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  SeqILDL<MatrixType,VectorType,VectorType> ildl(stiffnessMatrix,1.0);

  // Preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,
                          ildl, // preconditioner
                          1e-4, // desired residual reduction factor
                          100,   // maximum number of iterations
                          2);   // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  auto xFunction = Functions::makeDiscreteGlobalBasisFunction<double>(feBasis, x);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView, Dune::refinementLevels(2));
  vtkWriter.addVertexData(xFunction, VTK::FieldInfo("x", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("poisson-pq2");

 }
// Error handling
 catch (Exception& e) {
    std::cout << e.what() << std::endl;
 }
