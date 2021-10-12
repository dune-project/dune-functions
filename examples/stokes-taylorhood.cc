// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <array>
#include <vector>

#include <dune/common/bitsetvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/transpose.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#define BLOCKEDBASIS 1

// { using_namespace_dune_begin }
using namespace Dune;
// { using_namespace_dune_end }

// Compute the stiffness matrix for a single element
// { local_assembler_signature_begin }
template <class LocalView>
void getLocalMatrix(
        const LocalView& localView,
        Matrix<FieldMatrix<double,1,1>>& elementMatrix)
// { local_assembler_signature_end }
{
  // Get the grid element from the local FE basis view
  // { local_assembler_get_element_information_begin }
  using Element = typename LocalView::Element;
  const Element element = localView.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();
  // { local_assembler_get_element_information_end }

  // Set all matrix entries to zero
  // { initialize_element_matrix_begin }
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0;      // fills the entire matrix with zeros
  // { initialize_element_matrix_end }

  // Get set of shape functions for this element
  // { get_local_fe_begin }
  using namespace Indices;
  const auto& velocityLocalFiniteElement                   /*@\label{li:stokes_taylorhood_get_velocity_lfe}@*/
          = localView.tree().child(_0,0).finiteElement();
  const auto& pressureLocalFiniteElement
          = localView.tree().child(_1).finiteElement();    /*@\label{li:stokes_taylorhood_get_pressure_lfe}@*/
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
    // The transposed inverse Jacobian of the map from the
    // reference element to the element
    const auto jacobianInverseTransposed
            = geometry.jacobianInverseTransposed(quadPoint.position());

    // The multiplicative factor in the integral transformation formula
    const auto integrationElement
            = geometry.integrationElement(quadPoint.position());
    // { quad_loop_preamble_end }

    ///////////////////////////////////////////////////////////////////////
    //  Velocity--velocity coupling
    ///////////////////////////////////////////////////////////////////////

    // The gradients of the shape functions on the reference element
    // { velocity_gradients_begin }
    std::vector<FieldMatrix<double,1,dim> > referenceJacobians;
    velocityLocalFiniteElement.localBasis().evaluateJacobian(
            quadPoint.position(),
            referenceJacobians);

    // Compute the shape function gradients on the grid element
    std::vector<FieldMatrix<double,1,dim> > jacobians(referenceJacobians.size());
    for (size_t i=0; i<jacobians.size(); i++)
      jacobians[i] = referenceJacobians[i] * transpose(jacobianInverseTransposed);
    // { velocity_gradients_end }

    // Compute the actual matrix entries
    // { velocity_velocity_coupling_begin }
    for (size_t i=0; i<velocityLocalFiniteElement.size(); i++)
      for (size_t j=0; j<velocityLocalFiniteElement.size(); j++ )
        for (size_t k=0; k<dim; k++)
        {
          size_t row = localView.tree().child(_0,k).localIndex(i);                    /*@\label{li:stokes_taylorhood_compute_vv_element_matrix_row}@*/
          size_t col = localView.tree().child(_0,k).localIndex(j);                    /*@\label{li:stokes_taylorhood_compute_vv_element_matrix_column}@*/
          elementMatrix[row][col] += (jacobians[i] * transpose(jacobians[j]))
                                     * quadPoint.weight() * integrationElement;  /*@\label{li:stokes_taylorhood_update_vv_element_matrix}@*/
        }
    // { velocity_velocity_coupling_end }

    ///////////////////////////////////////////////////////////////////////
    //  Velocity--pressure coupling
    ///////////////////////////////////////////////////////////////////////

    // The values of the pressure shape functions
    // { pressure_values_begin }
    std::vector<FieldVector<double,1> > pressureValues;
    pressureLocalFiniteElement.localBasis().evaluateFunction(
            quadPoint.position(),
            pressureValues);
    // { pressure_values_end }

    // Compute the actual matrix entries
    // { velocity_pressure_coupling_begin }
    for (size_t i=0; i<velocityLocalFiniteElement.size(); i++)
      for (size_t j=0; j<pressureLocalFiniteElement.size(); j++ )
        for (size_t k=0; k<dim; k++)
        {
          size_t vIndex = localView.tree().child(_0,k).localIndex(i); /*@\label{li:stokes_taylorhood_compute_vp_element_matrix_row}@*/
          size_t pIndex = localView.tree().child(_1).localIndex(j);   /*@\label{li:stokes_taylorhood_compute_vp_element_matrix_column}@*/

          elementMatrix[vIndex][pIndex] +=                    /*@\label{li:stokes_taylorhood_update_vp_element_matrix_a}@*/
                  jacobians[i][0][k] * pressureValues[j]
                  * quadPoint.weight() * integrationElement;
          elementMatrix[pIndex][vIndex] +=
                  jacobians[i][0][k] * pressureValues[j]
                  * quadPoint.weight() * integrationElement;  /*@\label{li:stokes_taylorhood_update_vp_element_matrix_b}@*/
        }
    // { velocity_pressure_coupling_end }

  }

}


// Set the occupation pattern of the stiffness matrix
template <class Basis, class MatrixType>
void setOccupationPattern(const Basis& basis, MatrixType& matrix)
{
  enum {dim = Basis::GridView::dimension};

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  std::array<std::array<MatrixIndexSet, 2>, 2> nb;

  // Set sizes of the 2x2 submatrices
  for (size_t i=0; i<2; i++)
    for (size_t j=0; j<2; j++)
      nb[i][j].resize(basis.size({i}), basis.size({j}));

  // A view on the FE basis on a single element
  auto localView = basis.localView();

  // Loop over all leaf elements
  for(const auto& element : elements(basis.gridView()))
  {
    // Bind the local  view to the current element
    localView.bind(element);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<localView.size(); i++) {

      // Global index of the i-th local degree of freedom of the current element
      auto row = localView.index(i);

      for (size_t j=0; j<localView.size(); j++ ) {

        // Global index of the j-th local degree of freedom of the current element
        auto col = localView.index(j);

        nb[row[0]][col[0]].add(row[1],col[1]);

      }

    }

  }

  // Give the matrix the occupation pattern we want.
  using namespace Indices;
#if !BLOCKEDBASIS
  matrix.setSize(2,2);
#endif
  nb[0][0].exportIdx(matrix[_0][_0]);
  nb[0][1].exportIdx(matrix[_0][_1]);
  nb[1][0].exportIdx(matrix[_1][_0]);
  nb[1][1].exportIdx(matrix[_1][_1]);
}


#if BLOCKEDBASIS
// { matrixentry_begin }
template<class Matrix, class MultiIndex>
decltype(auto) matrixEntry(
        Matrix& matrix, const MultiIndex& row, const MultiIndex& col)
{
  using namespace Indices;
  if ((row[0]==0) and (col[0]==0))
    return matrix[_0][_0][row[1]][col[1]][row[2]][col[2]];
  if ((row[0]==0) and (col[0]==1))
    return matrix[_0][_1][row[1]][col[1]][row[2]][0];
  if ((row[0]==1) and (col[0]==0))
    return matrix[_1][_0][row[1]][col[1]][0][col[2]];
  return matrix[_1][_1][row[1]][col[1]][0][0];  /*@\label{li:matrixentry_pressure_pressure}@*/
}
// { matrixentry_end }
#else
template<class Matrix, class MultiIndex>
decltype(auto) matrixEntry(Matrix& matrix, const MultiIndex& row, const MultiIndex& col)
{
  return matrix[row[0]][col[0]][row[1]][col[1]];
}
#endif


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
// { global_assembler_signature_begin }
template <class Basis, class MatrixType>
void assembleStokesMatrix(const Basis& basis, MatrixType& matrix)
// { global_assembler_signature_end }
{
  // { setup_matrix_pattern_begin }
  // Set matrix size and occupation pattern
  setOccupationPattern(basis, matrix);

  // Set all entries to zero
  matrix = 0;
  // { setup_matrix_pattern_end }

  // A view on the FE basis on a single element
  // { get_localview_begin }
  auto localView     = basis.localView();
  // { get_localview_end }

  // A loop over all elements of the grid
  // { element_loop_and_bind_begin }
  for (const auto& element : elements(basis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    // { element_loop_and_bind_end }

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    // { setup_element_stiffness_begin }
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix);
    // { setup_element_stiffness_end }

    // Add element stiffness matrix onto the global stiffness matrix
    // { accumulate_global_matrix_begin }
    for (size_t i=0; i<elementMatrix.N(); i++)
    {
      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localView.index(i);                /*@\label{li:stokes_taylorhood_get_global_row_index}@*/

      for (size_t j=0; j<elementMatrix.M(); j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localView.index(j);                /*@\label{li:stokes_taylorhood_get_global_column_index}@*/
        matrixEntry(matrix, row, col) += elementMatrix[i][j];  /*@\label{li:stokes_taylorhood_scatter_matrix_indices}@*/
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
  // { mpi_setup_end }

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  // { grid_setup_begin }
  const int dim = 2;
  using GridType = YaspGrid<dim>;
  FieldVector<double,dim> upperRight = {1, 1};
  std::array<int,dim> elements = {{4, 4}};
  GridType grid(upperRight,elements);

  using GridView = typename GridType::LeafGridView;
  GridView gridView = grid.leafGridView();
  // { grid_setup_end }

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

#if BLOCKEDBASIS
  // { function_space_basis_begin }
  using namespace Functions::BasisFactory;

  constexpr std::size_t p = 1; // pressure order for Taylor-Hood

  auto taylorHoodBasis = makeBasis(
          gridView,
          composite(
            power<dim>(
              lagrange<p+1>(),
              blockedInterleaved()),
            lagrange<p>()
          ));
  // { function_space_basis_end }
#else
  using namespace Functions::BasisFactory;

  static const std::size_t p = 1; // pressure order for Taylor-Hood
  auto taylorHoodBasis = makeBasis(
          gridView,
          composite(
            power<dim>(
              lagrange<p+1>(),
              flatInterleaved()),
            lagrange<p>()
          ));
#endif



  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

#if BLOCKEDBASIS
  // { linear_algebra_setup_begin }
  using VelocityVector = BlockVector<FieldVector<double,dim>>;
  using PressureVector = BlockVector<FieldVector<double,1>>;
  using VectorType = MultiTypeBlockVector<VelocityVector, PressureVector>;

  using VelocityBitVector = BlockVector<FieldVector<char,dim>>;
  using PressureBitVector = BlockVector<FieldVector<char,1>>;
  using BitVectorType = MultiTypeBlockVector<VelocityBitVector, PressureBitVector>;

  using Matrix00 = BCRSMatrix<FieldMatrix<double,dim,dim>>;
  using Matrix01 = BCRSMatrix<FieldMatrix<double,dim,1>>;
  using Matrix10 = BCRSMatrix<FieldMatrix<double,1,dim>>;
  using Matrix11 = BCRSMatrix<FieldMatrix<double,1,1>>;   /*@\label{li:matrix_type_pressure_pressure}@*/
  using MatrixRow0 = MultiTypeBlockVector<Matrix00, Matrix01>;
  using MatrixRow1 = MultiTypeBlockVector<Matrix10, Matrix11>;
  using MatrixType = MultiTypeBlockMatrix<MatrixRow0,MatrixRow1>;
  // { linear_algebra_setup_end }
#else
  using VectorType = BlockVector<BlockVector<FieldVector<double,1> > >;
  using BitVectorType = BlockVector<BlockVector<FieldVector<char,1> > >;
  using MatrixType = Matrix<BCRSMatrix<FieldMatrix<double,1,1> > >;
#endif

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  // { rhs_assembly_begin }
  VectorType rhs;

  auto rhsBackend = Dune::Functions::istlVectorBackend(rhs);

  rhsBackend.resize(taylorHoodBasis);
  rhs = 0;                                 /*@\label{li:stokes_taylorhood_set_rhs_to_zero}@*/
  // { rhs_assembly_end }

  // { matrix_assembly_begin }
  MatrixType stiffnessMatrix;
  assembleStokesMatrix(taylorHoodBasis, stiffnessMatrix);   /*@\label{li:stokes_taylorhood_call_to_assemblestokesmatrix}@*/
  // { matrix_assembly_end }

  /////////////////////////////////////////////////////////
  // Set Dirichlet values.
  // Only velocity components have Dirichlet boundary values
  /////////////////////////////////////////////////////////

  // { initialize_boundary_dofs_vector_begin }

  BitVectorType isBoundary;

  auto isBoundaryBackend = Dune::Functions::istlVectorBackend(isBoundary);
  isBoundaryBackend.resize(taylorHoodBasis);
  isBoundary = false;
  // { initialize_boundary_dofs_vector_end }

  // { determine_boundary_dofs_begin }
  using namespace Indices;
  Functions::forEachBoundaryDOF(
          Functions::subspaceBasis(taylorHoodBasis, _0),
          [&] (auto&& index) {
            isBoundaryBackend[index] = true;
          });
  // { determine_boundary_dofs_end }

  // { interpolate_dirichlet_values_begin }
  using Coordinate = GridView::Codim<0> ::Geometry::GlobalCoordinate;
  using VelocityRange = FieldVector<double,dim>;
  auto&& velocityDirichletData = [](Coordinate x)
  {
    return VelocityRange{0.0, double(x[0] < 1e-8)};
  };

  Functions::interpolate(
          Functions::subspaceBasis(taylorHoodBasis, _0), rhs,
          velocityDirichletData,
          isBoundary);
  // { interpolate_dirichlet_values_end }

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  // { set_dirichlet_matrix_begin }
  auto localView = taylorHoodBasis.localView();
  for(const auto& element : Dune::elements(taylorHoodBasis.gridView()))
  {
    localView.bind(element);
    for (size_t i=0; i<localView.size(); ++i)
    {
      auto row = localView.index(i);
      // If row corresponds to a boundary entry, modify
      // it to be an identity matrix row
      if (isBoundaryBackend[row])
        for (size_t j=0; j<localView.size(); ++j)
        {
          auto col = localView.index(j);
          matrixEntry(stiffnessMatrix, row, col) = (i==j) ? 1 : 0;
        }
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
  MatrixAdapter<MatrixType,VectorType,VectorType> stiffnessOperator(stiffnessMatrix);

  // Fancy (but only) way to not have a preconditioner at all
  Richardson<VectorType,VectorType> preconditioner(1.0);

  // Construct the actual iterative solver
  RestartedGMResSolver<VectorType> solver(
          stiffnessOperator,  // operator to invert
          preconditioner,     // preconditioner for interation
          1e-10,              // desired residual reduction factor
          500,                // number of iterations between restarts
          500,                // maximum number of iterations
          2);                 // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  solver.apply(x, rhs, statistics);
  // { stokes_solve_end }

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  // { make_result_functions_begin }
  using VelocityRange = FieldVector<double,dim>;
  using PressureRange = double;

  auto velocityFunction
          = Functions::makeDiscreteGlobalBasisFunction<VelocityRange>(
            Functions::subspaceBasis(taylorHoodBasis, _0), x);
  auto pressureFunction
          = Functions::makeDiscreteGlobalBasisFunction<PressureRange>(
            Functions::subspaceBasis(taylorHoodBasis, _1), x);
  // { make_result_functions_end }

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  // { vtk_output_begin }
  SubsamplingVTKWriter<GridView> vtkWriter(
          gridView,
          refinementLevels(2));
  vtkWriter.addVertexData(
          velocityFunction,
          VTK::FieldInfo("velocity", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.addVertexData(
          pressureFunction,
          VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("stokes-taylorhood-result");
  // { vtk_output_end }

 }
// Error handling
 catch (Exception& e) {
    std::cout << e.what() << std::endl;
 }
