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
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

// { using_namespace_dune_begin }
using namespace Dune;
// { using_namespace_dune_end }

// Compute the stiffness matrix for a single element
// { local_assembler_signature_begin }
template <class LocalView>
void getLocalMatrix(const LocalView& localView,
                    Matrix<FieldMatrix<double,1,1> >& elementMatrix)
// { local_assembler_signature_end }
{
  // Get the grid element from the local FE basis view
  // { local_assembler_get_element_information_begin }
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();
  // { local_assembler_get_element_information_end }

  // Set all matrix entries to zero
  // { initialize_element_matrix_begin }
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0;      // fills the entire matrix with zeroes
  // { initialize_element_matrix_end }

  // Get set of shape functions for this element
  // { get_local_fe_begin }
  using namespace Dune::TypeTree::Indices;
  const auto& velocityLocalFiniteElement = localView.tree().child(_0,0).finiteElement();  /*@\label{li:stokes_taylorhood_get_velocity_lfe}@*/
  const auto& pressureLocalFiniteElement = localView.tree().child(_1).finiteElement();    /*@\label{li:stokes_taylorhood_get_pressure_lfe}@*/
  // { get_local_fe_end }

  // Get a quadrature rule
  // { begin_quad_loop_begin }
  int order = 2*(dim*velocityLocalFiniteElement.localBasis().order()-1);
  const auto& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad)
  {
    // { begin_quad_loop_end }
    // { quad_loop_preamble_begin }
    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto jacobian = geometry.jacobianInverseTransposed(quadPoint.position());

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPoint.position());
    // { quad_loop_preamble_end }

    ///////////////////////////////////////////////////////////////////////
    //  Velocity--velocity coupling
    ///////////////////////////////////////////////////////////////////////

    // The gradients of the shape functions on the reference element
    // { velocity_gradients_begin }
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    velocityLocalFiniteElement.localBasis().evaluateJacobian(quadPoint.position(), referenceGradients);

    // Compute the shape function gradients on the grid element
    std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
    for (size_t i=0; i<gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);
    // { velocity_gradients_end }

    // Compute the actual matrix entries
    // { velocity_velocity_coupling_begin }
    for (size_t i=0; i<velocityLocalFiniteElement.size(); i++)
      for (size_t j=0; j<velocityLocalFiniteElement.size(); j++ )
        for (size_t k=0; k<dim; k++)
        {
          size_t row = localView.tree().child(_0,k).localIndex(i);                    /*@\label{li:stokes_taylorhood_compute_vv_element_matrix_row}@*/
          size_t col = localView.tree().child(_0,k).localIndex(j);                    /*@\label{li:stokes_taylorhood_compute_vv_element_matrix_column}@*/
          elementMatrix[row][col] += ( gradients[i] * gradients[j] ) * quadPoint.weight() * integrationElement;  /*@\label{li:stokes_taylorhood_update_vv_element_matrix}@*/
        }
    // { velocity_velocity_coupling_end }

    ///////////////////////////////////////////////////////////////////////
    //  Velocity--pressure coupling
    ///////////////////////////////////////////////////////////////////////

    // The values of the pressure shape functions
    // { pressure_values_begin }
    std::vector<FieldVector<double,1> > pressureValues;
    pressureLocalFiniteElement.localBasis().evaluateFunction(quadPoint.position(), pressureValues);
    // { pressure_values_end }

    // Compute the actual matrix entries
    // { velocity_pressure_coupling_begin }
    for (size_t i=0; i<velocityLocalFiniteElement.size(); i++)
      for (size_t j=0; j<pressureLocalFiniteElement.size(); j++ )
        for (size_t k=0; k<dim; k++)
        {
          size_t vIndex = localView.tree().child(_0,k).localIndex(i);                      /*@\label{li:stokes_taylorhood_compute_vp_element_matrix_row}@*/
          size_t pIndex = localView.tree().child(_1).localIndex(j);                        /*@\label{li:stokes_taylorhood_compute_vp_element_matrix_column}@*/

          elementMatrix[vIndex][pIndex] += gradients[i][k] * pressureValues[j] * quadPoint.weight() * integrationElement;  /*@\label{li:stokes_taylorhood_update_vp_element_matrix_a}@*/
          elementMatrix[pIndex][vIndex] += gradients[i][k] * pressureValues[j] * quadPoint.weight() * integrationElement;  /*@\label{li:stokes_taylorhood_update_vp_element_matrix_b}@*/
        }
    // { velocity_pressure_coupling_end }

  }

}


// Get the occupation pattern of the stiffness matrix
template <class Basis>
void getOccupationPattern(const Basis& basis,
                          std::array<std::array<MatrixIndexSet,2>,2>& nb)
{
  enum {dim = Basis::GridView::dimension};

  // Set sizes of the 2x2 submatrices
  for (size_t i=0; i<2; i++)
    for (size_t j=0; j<2; j++)
      nb[i][j].resize(basis.size({i}), basis.size({j}));

  // A view on the FE basis on a single element
  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();

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
// { global_assembler_signature_begin }
template <class Basis>
void assembleStokesMatrix(const Basis& basis,
                          Matrix<BCRSMatrix<FieldMatrix<double,1,1> > >& matrix)
// { global_assembler_signature_end }
{
  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  // { setup_matrix_pattern_begin }
  array<array<MatrixIndexSet, 2>, 2> occupationPattern;
  getOccupationPattern(basis, occupationPattern);

  // ... and give it the occupation pattern we want.
  matrix.setSize(2,2);
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      occupationPattern[i][j].exportIdx(matrix[i][j]);  /*@\label{li:stokes_taylorhood_setup_matrix_patterns}@*/

  // Set all entries to zero
  matrix = 0;                           /*@\label{li:stokes_taylorhood_set_matrix_to_zero}@*/
  // { setup_matrix_pattern_end }

  // A view on the FE basis on a single element
  // { get_localview_begin }
  auto localView     = basis.localView();
  auto localIndexSet = basis.localIndexSet();
  // { get_localview_end }

  // A loop over all elements of the grid
  // { element_loop_and_bind_begin }
  for (const auto& element : elements(basis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localIndexSet.bind(localView);
    // { element_loop_and_bind_end }

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    // { setup_element_stiffness_begin }
    Matrix<FieldMatrix<double,1,1> > elementMatrix;       /*@\label{li:stokes_taylorhood_select_element_matrix_type}@*/
    getLocalMatrix(localView, elementMatrix);
    // { setup_element_stiffness_end }

    // Add element stiffness matrix onto the global stiffness matrix
    // { accumulate_global_matrix_begin }
    for (size_t i=0; i<elementMatrix.N(); i++)
    {
      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localIndexSet.index(i);                /*@\label{li:stokes_taylorhood_get_global_row_index}@*/

      for (size_t j=0; j<elementMatrix.M(); j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localIndexSet.index(j);                /*@\label{li:stokes_taylorhood_get_global_column_index}@*/
        matrix[row[0]][col[0]][row[1]][col[1]] += elementMatrix[i][j];  /*@\label{li:stokes_taylorhood_scatter_matrix_indices}@*/
      }

    }
    // { accumulate_global_matrix_end }

  }

}



// { main_begin }
int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> bbox = {1, 1};
  std::array<int,dim> elements = {4, 4};
  GridType grid(bbox,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();
  // { grid_setup_end }

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  // { function_space_basis_begin }
//  typedef Functions::TaylorHoodBasis<GridView> TaylorHoodBasis;        /*@\label{li:stokes_taylorhood_select_taylorhoodbasis}@*/
//  TaylorHoodBasis taylorHoodBasis(gridView);

  using namespace Functions::BasisBuilder;

  static const std::size_t K = 1; // pressure order for Taylor-Hood

  auto taylorHoodBasis = makeBasis(
    gridView,
    composite(
      power<dim>(
        lagrange<K+1>(),
        flatInterleaved()),
      lagrange<K>()
    ));


  // { function_space_basis_end }

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  // { linear_algebra_setup_begin }
  typedef BlockVector<BlockVector<FieldVector<double,1> > > VectorType;
  typedef Matrix<BCRSMatrix<FieldMatrix<double,1,1> > > MatrixType;
  // { linear_algebra_setup_end }

  // { rhs_assembly_begin }
  VectorType rhs;

  typedef Dune::Functions::HierarchicVectorWrapper<VectorType, double> HierarchicVectorView;
  HierarchicVectorView(rhs).resize(taylorHoodBasis);

  rhs = 0;                                 /*@\label{li:stokes_taylorhood_set_rhs_to_zero}@*/
  // { rhs_assembly_end }

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////
  // { matrix_assembly_begin }
  MatrixType stiffnessMatrix;
  assembleStokesMatrix(taylorHoodBasis, stiffnessMatrix);   /*@\label{li:stokes_taylorhood_call_to_assemblestokesmatrix}@*/
  // { matrix_assembly_end }

  // Set Dirichlet values
  // Only velocity components have Dirichlet boundary values
  // { boundary_predicate_begin }
  using Coordinate = GridView::Codim<0> ::Geometry::GlobalCoordinate;

  auto boundaryIndicator = [&bbox](Coordinate x) {
    bool isBoundary = false;
    for (std::size_t j=0; j<x.size(); j++)
      isBoundary |= x[j] < 1e-8 || x[j] > bbox[j] - 1e-8;
    return isBoundary;
  };
  // { boundary_predicate_end }

  // { interpolate_boundary_predicate_begin }
  using namespace TypeTree::Indices;

  using BitVectorType = BlockVector<BlockVector<FieldVector<char,1> > >;
  using HierarchicBitVectorView = Functions::HierarchicVectorWrapper<BitVectorType, char>;

  BitVectorType isBoundary;

  for(int i=0; i<dim; ++i)
  {
    auto velocityComponentSpace = Functions::subspaceBasis(taylorHoodBasis, _0, i);
    interpolate(velocityComponentSpace, HierarchicBitVectorView(isBoundary), boundaryIndicator);
  }
  // { interpolate_boundary_predicate_end }

  // { interpolate_dirichlet_values_begin }
  typedef FieldVector<double,dim> VelocityRange;
  auto&& velocityDirichletData = [](Coordinate x)->VelocityRange { return {0.0, double(x[0] < 1e-8)};};

  interpolate(Functions::subspaceBasis(taylorHoodBasis, _0),
                        HierarchicVectorView(rhs),
                        velocityDirichletData,
                        HierarchicBitVectorView(isBoundary));
  // { interpolate_dirichlet_values_end }

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  // { set_dirichlet_matrix_begin }
  for (size_t i=0; i<stiffnessMatrix[0][0].N(); i++)
  {
    if (isBoundary[0][i])
    {
      // Upper left matrix block
      auto cIt    = stiffnessMatrix[0][0][i].begin();
      auto cEndIt = stiffnessMatrix[0][0][i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? 1 : 0;

      // Upper right matrix block
      for(auto&& entry: stiffnessMatrix[0][1][i])
        entry = 0.0;
    }

  }
  // { set_dirichlet_matrix_end }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////
  // { stokes_solve_begin }
  // Start from the rhs vector; that way the Dirichlet entries are already correct
  VectorType x = rhs;

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Fancy (but only) way to not have a preconditioner at all
  Richardson<VectorType,VectorType> preconditioner(1.0);

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
  // { stokes_solve_end }

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  // { stokes_output_begin }
  using VelocityRange = FieldVector<double,dim>;
  using PressureRange = double;

  auto velocityFunction = Functions::makeDiscreteGlobalBasisFunction<VelocityRange>(Functions::subspaceBasis(taylorHoodBasis, _0),
                                                                                    HierarchicVectorView(x));
  auto pressureFunction = Functions::makeDiscreteGlobalBasisFunction<PressureRange>(Functions::subspaceBasis(taylorHoodBasis, _1),
                                                                                    HierarchicVectorView(x));

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(velocityFunction, VTK::FieldInfo("velocity", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.addVertexData(pressureFunction, VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("stokes-taylorhood-result");
  // { stokes_output_end }

 }
// Error handling
 catch (Exception e) {
    std::cout << e << std::endl;
 }
