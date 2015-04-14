// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/functionspacebases/pq2nodalbasis.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

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

  // Get set of shape functions for this element
  const auto& velocityLocalFiniteElement = localView.tree().template child<0>().child(0).finiteElement();
  const auto& pressureLocalFiniteElement = localView.tree().template child<1>().finiteElement();

  // Set all matrix entries to zero
  elementMatrix.setSize(dim*velocityLocalFiniteElement.size() + pressureLocalFiniteElement.size(),
                        dim*velocityLocalFiniteElement.size() + pressureLocalFiniteElement.size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // Get a quadrature rule
  int order = 2*(dim*velocityLocalFiniteElement.localBasis().order()-1);
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (size_t pt=0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    ///////////////////////////////////////////////////////////////////////
    //  Velocity--velocity coupling
    ///////////////////////////////////////////////////////////////////////

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    velocityLocalFiniteElement.localBasis().evaluateJacobian(quadPos, referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
    for (size_t i=0; i<gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);

    // Compute the actual matrix entries
    for (size_t i=0; i<velocityLocalFiniteElement.size(); i++)
      for (size_t j=0; j<velocityLocalFiniteElement.size(); j++ )
        for (size_t k=0; k<dim; k++)
          elementMatrix[dim*i+k][dim*j+k] += ( gradients[i] * gradients[j] ) * quad[pt].weight() * integrationElement;

    ///////////////////////////////////////////////////////////////////////
    //  Velocity--pressure coupling
    ///////////////////////////////////////////////////////////////////////

    // The values of the pressure shape functions
    std::vector<FieldVector<double,1> > pressureValues;
    pressureLocalFiniteElement.localBasis().evaluateFunction(quadPos, pressureValues);

    // Compute the actual matrix entries
    for (size_t i=0; i<velocityLocalFiniteElement.size(); i++)
      for (size_t j=0; j<pressureLocalFiniteElement.size(); j++ )
        for (size_t k=0; k<dim; k++)
        {
          elementMatrix[dim*i+k][dim*velocityLocalFiniteElement.size() + j] += gradients[i][k] * pressureValues[j] * quad[pt].weight() * integrationElement;
          elementMatrix[dim*velocityLocalFiniteElement.size() + j][dim*i+k] += gradients[i][k] * pressureValues[j] * quad[pt].weight() * integrationElement;
        }

  }

}


// Get the occupation pattern of the stiffness matrix
template <class Basis>
void getOccupationPattern(const Basis& basis,
                          std::array<std::array<MatrixIndexSet,2>,2>& nb)
{
  enum {dim = Basis::GridView::dimension};

  // Total number of grid vertices
  auto basisIndexSet = basis.indexSet();

  for (size_t i=0; i<2; i++)
    for (size_t j=0; j<2; j++)
      nb[i][j].resize(basisIndexSet.size({i}), basisIndexSet.size({j}));

  // A view on the FE basis on a single element
  typename Basis::LocalView localView(&basis);
  auto localIndexSet = basisIndexSet.localIndexSet();

  // Loop over all leaf elements
  for(const auto& e : elements(basis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(e);
    localIndexSet.bind(localView);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<localIndexSet.size(); i++) {

      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<localIndexSet.size(); j++ ) {

        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localIndexSet.index(j);
        nb[row[0]][col[0]].add(row[1],col[1]);

      }

    }

  }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class Basis>
void assembleStokesProblem(const Basis& basis,
                           Matrix<BCRSMatrix<FieldMatrix<double,1,1> > >& matrix)
{
  // Get the grid view from the finite element basis
  typedef typename Basis::GridView GridView;
  GridView gridView = basis.gridView();

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  array<array<MatrixIndexSet, 2>, 2> occupationPattern;
  getOccupationPattern(basis, occupationPattern);

  // ... and give it the occupation pattern we want.
  matrix.setSize(2,2);
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      occupationPattern[i][j].exportIdx(matrix[i][j]);

  // Set all entries to zero
  matrix = 0;

  // A view on the FE basis on a single element
  auto localView     = basis.localView();
  auto localIndexSet = basis.indexSet().localIndexSet();

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(e);
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


// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
template <class VelocityBasis>
void boundaryTreatment (const VelocityBasis& velocityBasis,
                        const FieldVector<double,VelocityBasis::GridView::dimension>& bbox,
                        std::vector<int>& dirichletNodes)
{
  static const int dim = VelocityBasis::GridView::dimension;

  // Interpolating the identity function wrt to a Lagrange basis
  // yields the positions of the Lagrange nodes
  BlockVector<FieldVector<double,dim> > lagrangeNodes;
  interpolate(velocityBasis, lagrangeNodes, [](FieldVector<double,dim> x){ return x; });

  dirichletNodes.resize(dim*velocityBasis.indexSet().size());

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    for (int j=0; j<dim; j++)
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8 || lagrangeNodes[i][j] > bbox[j]-1e-8;

    if (isBoundary)
    {
      for (int j=0; j<dim; j++)
        dirichletNodes[dim*i+j] = 1;

      if (lagrangeNodes[i][0] < 1e-8)
        dirichletNodes[dim*i+1] = 2;
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
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {4, 4};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::TaylorHoodBasis<GridView> TaylorHoodBasis;
  TaylorHoodBasis taylorHoodBasis(gridView);

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  typedef BlockVector<BlockVector<FieldVector<double,1> > > VectorType;
  typedef Matrix<BCRSMatrix<FieldMatrix<double,1,1> > > MatrixType;

  VectorType rhs(2);
  rhs[0].resize(taylorHoodBasis.indexSet().size({0}));
  rhs[1].resize(taylorHoodBasis.indexSet().size({1}));
  rhs[0] = 0;
  rhs[1] = 0;
  MatrixType stiffnessMatrix;

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  assembleStokesProblem(taylorHoodBasis, stiffnessMatrix);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(2);
  x[0].resize(taylorHoodBasis.indexSet().size({0}));
  x[1].resize(taylorHoodBasis.indexSet().size({1}));
  x[0] = 0;
  x[1] = 0;

  // Determine Dirichlet dofs
  typedef Functions::PQ2NodalBasis<GridView> VelocityBasis;
  VelocityBasis velocityBasis(gridView);

  typedef Functions::PQ1NodalBasis<GridView> PressureBasis;
  PressureBasis pressureBasis(gridView);

  std::vector<int> dirichletNodes;
  boundaryTreatment(velocityBasis, l, dirichletNodes);

  // Set Dirichlet values
  // Only velocity components have Dirichlet boundary values
  for (size_t i=0; i<rhs[0].size(); i++)
  {
    if (dirichletNodes[i]==1)
      // The zero is the value of the Dirichlet boundary condition
      rhs[0][i] = 0;
    else if (dirichletNodes[i]==2)
      rhs[0][i] = 1;
  }

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<stiffnessMatrix[0][0].N(); i++) {

    if (dirichletNodes[i]) {

      // Upper left matrix block
      auto cIt    = stiffnessMatrix[0][0][i].begin();
      auto cEndIt = stiffnessMatrix[0][0][i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? 1 : 0;

      // Upper right matrix block
      stiffnessMatrix[0][1][i] = 0;

    }

  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Fancy (but only) way to not have a preconditioner at all
  Richardson<VectorType,VectorType> preconditioner(1.0);

  // Preconditioned conjugate-gradient solver
  RestartedGMResSolver<VectorType> solver(op,
                                          preconditioner,
                                          1e-4,  // desired residual reduction factor
                                          5,     // number of iterations between restarts
                                          100,   // maximum number of iterations
                                          2);    // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  solver.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  typedef BlockVector<FieldVector<double,dim> > VelocityVectorType;
  typedef BlockVector<FieldVector<double,1> >   PressureVectorType;

  VelocityVectorType velocity(velocityBasis.indexSet().size());
  for (size_t i=0; i<velocity.size(); i++)
    for (int j=0; j<dim; j++)
      velocity[i][j] = x[0][dim*i+j];

  PressureVectorType pressure(pressureBasis.indexSet().size());
  for (size_t i=0; i<pressure.size(); i++)
    pressure[i] = x[1][i];


  Dune::Functions::DiscreteScalarGlobalBasisFunction<VelocityBasis,VelocityVectorType> velocityFunction(velocityBasis,velocity);
  auto localVelocityFunction = localFunction(velocityFunction);

  Dune::Functions::DiscreteScalarGlobalBasisFunction<PressureBasis,PressureVectorType> pressureFunction(pressureBasis,pressure);
  auto localPressureFunction = localFunction(pressureFunction);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(localVelocityFunction, VTK::FieldInfo("velocity", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.addVertexData(localPressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("functions-stokes");

 }
// Error handling
 catch (Exception e) {
    std::cout << e << std::endl;
 }
