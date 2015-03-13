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
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/pq2nodalbasis.hh>

using namespace Dune;

// HACK: THE FOLLOWING CLASS SHOULD REALLY BE AVAILABLE FROM THE DUNE-FUNCTIONS MODULE
//  However, currently it isn't, so I add it here to get a preliminary working program.


/** \brief A VTK basis grid function.
 *
 *  This function "evaluates" by evaluating the global basis and interpolating the function values.
 *  \tparam Basis   The global basis.
 *  \tparam CoefficientType The vector type for the coefficients.
*/
template<class Basis, class CoefficientType>
class VTKBasisGridFunction : public Dune::VTKFunction<typename Basis::GridView>
{
  typedef Dune::VTKFunction<typename Basis::GridView> Base;
  typedef typename  CoefficientType::value_type RangeType;
public:
  typedef typename Base::Entity Entity;
  typedef typename Base::ctype ctype;
  using Base::dim;

  /** \brief Construct from given global basis, coefficient vector and name.
   *
   *  \param basis    The global basis.
   *  \param v    A corresponding vector of coefficients.
   *  \param s    A name of the function.
   */
  VTKBasisGridFunction(const Basis &basis, const CoefficientType &v, const std::string &s) :
      basis_(basis),
      coeffs_(v),
      s_( s )
  {
    if (v.size() !=basis_.indexSet().size())
      DUNE_THROW(Dune::IOError, "VTKGridFunction: Coefficient vector is not matching the basis");
  }

  /** \brief Get the number of components the function has. */
  virtual int ncomps () const
  {
    return CoefficientType::value_type::dimension;
  }

  /** \brief Locally evaluate a component of the function.
   *
   *  \param comp The component to evaluate.
   *  \param e    The element the local coordinates are taken from.
   *  \param xi   The local coordinates where to evaluate the function.
   */
  virtual double evaluate (int comp, const Entity &e,
                           const Dune::FieldVector<ctype,dim> &xi) const
  {
    typename Basis::LocalView localView(&basis_);
    localView.bind(e);
    auto basisIndexSet = basis_.indexSet();
    auto localIndexSet = basisIndexSet.localIndexSet();
    localIndexSet.bind(localView);

    std::vector<RangeType> shapeFunctionValues;
    auto& basis = localView.tree().finiteElement().localBasis();
    basis.evaluateFunction(xi,shapeFunctionValues);
    RangeType r = 0;
    for (size_t i = 0; i < basis.size(); ++i)
      r += coeffs_[localIndexSet.index(i)[0]] * shapeFunctionValues[i];

    return r[comp];
  }

  /** \brief Get the name of that function. */
  virtual std::string name () const
  {
    return s_;
  }

  /** \brief Destructor. */
  virtual ~VTKBasisGridFunction() {}

private:
  const Basis &basis_;
  const CoefficientType &coeffs_;
  std::string s_;
};

// THIS IS THE END OF THE BLOCK OF TEMPORARY HELPER CODE.
// ALL CODE AFTER THIS LINE IS PART OF THE EXAMPLE PROPER.

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
template <class LocalView>
void getVolumeTerm( const LocalView& localView,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    const Dune::VirtualFunction<FieldVector<double,LocalView::Element::dimension>, double>* volumeTerm)
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

    double functionValue;
    volumeTerm->evaluate(element.geometry().global(quadPos), functionValue);

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
  auto basisIndexSet = feBasis.indexSet();
  auto n = basisIndexSet.size();

  nb.resize(n, n);

  // A view on the FE basis on a single element
  typename FEBasis::LocalView localView(&feBasis);
  auto localIndexSet = basisIndexSet.localIndexSet();

  // Loop over all leaf elements
  auto it    = feBasis.gridView().template begin<0>();
  auto endIt = feBasis.gridView().template end<0>  ();

  for (; it!=endIt; ++it)
  {
    // Bind the local FE basis view to the current element
    localView.bind(*it);
    localIndexSet.bind(localView);

    // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
    for (size_t i=0; i<localView.tree().size(); i++) {

      for (size_t j=0; j<localView.tree().size(); j++) {

        auto iIdx = localIndexSet.index(i)[0];
        auto jIdx = localIndexSet.index(j)[0];

        // Add a nonzero entry to the matrix
        nb.add(iIdx, jIdx);

      }

    }

  }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class FEBasis>
void assembleLaplaceMatrix(const FEBasis& feBasis,
                           BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                           BlockVector<FieldVector<double,1> >& rhs,
                           const VirtualFunction<Dune::FieldVector<double,FEBasis::GridView::dimension>,double>* volumeTerm)
{
  // Get the grid view from the finite element basis
  typedef typename FEBasis::GridView GridView;
  GridView gridView = feBasis.gridView();

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  getOccupationPattern(feBasis, occupationPattern);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  auto basisIndexSet = feBasis.indexSet();
  rhs.resize(basisIndexSet.size());

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // A view on the FE basis on a single element
  typename FEBasis::LocalView localView(&feBasis);
  auto localIndexSet = basisIndexSet.localIndexSet();

  // A loop over all elements of the grid
  auto it    = gridView.template begin<0>();
  auto endIt = gridView.template end<0>  ();

  for( ; it != endIt; ++it ) {

    // Bind the local FE basis view to the current element
    localView.bind(*it);
    localIndexSet.bind(localView);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++) {

      // The global index of the i-th local degree of freedom of the element 'it'
      auto row = localIndexSet.index(i)[0];

      for (size_t j=0; j<elementMatrix.M(); j++ ) {

        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localIndexSet.index(j)[0];
        matrix[row][col] += elementMatrix[i][j];

      }

    }

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    getVolumeTerm(localView, localRhs, volumeTerm);

    for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'it'
      auto row = localIndexSet.index(i)[0];
      rhs[row] += localRhs[i];

    }

  }

}


// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
template <class FEBasis>
void boundaryTreatment (const FEBasis& feBasis,
                        const FieldVector<double,FEBasis::GridView::dimension>& bbox,
                        std::vector<bool>& dirichletNodes )
{
  static const int dim = FEBasis::GridView::dimension;

  // Interpolating the identity function wrt to a Lagrange basis
  // yields the positions of the Lagrange nodes

  // TODO: We are hacking our way around the fact that interpolation
  // of vector-value functions is not supported yet.
  BlockVector<FieldVector<double,1> > lagrangeNodes0;
  BlockVector<FieldVector<double,1> > lagrangeNodes1;

  interpolate(feBasis, lagrangeNodes0, [](FieldVector<double,dim> x){ return x[0]; });
  interpolate(feBasis, lagrangeNodes1, [](FieldVector<double,dim> x){ return x[1]; });

  BlockVector<FieldVector<double,dim> > lagrangeNodes(lagrangeNodes0.size());
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    lagrangeNodes[i][0] = lagrangeNodes0[i];
    lagrangeNodes[i][1] = lagrangeNodes1[i];
  }

  dirichletNodes.resize(lagrangeNodes.size());

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    for (int j=0; j<dim; j++)
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8 || lagrangeNodes[i][j] > bbox[j]-1e-8;

    if (isBoundary)
      dirichletNodes[i] = true;
  }
}


// A class implementing the analytical right hand side.  Here simply constant '1'
template <int dim>
class RightHandSide
    : public VirtualFunction<FieldVector<double,dim>, double >
{
public:
    void evaluate(const FieldVector<double,dim>& in, double& out) const {
        out = 1;
    }
};



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
  std::array<int,dim> elements = {10, 10};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::PQ2NodalBasis<GridView> FEBasis;
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

  RightHandSide<dim> rightHandSide;
  assembleLaplaceMatrix(feBasis, stiffnessMatrix, rhs, &rightHandSide);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasis.indexSet().size());
  x = 0;

  // Determine Dirichlet dofs
  std::vector<bool> dirichletNodes;
  boundaryTreatment(feBasis, l, dirichletNodes);

  // Set Dirichlet values
  for (size_t i=0; i<rhs.size(); i++)
    if (dirichletNodes[i])
      // The zero is the value of the Dirichlet boundary condition
      rhs[i] = 0;

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<stiffnessMatrix.N(); i++) {

    if (dirichletNodes[i]) {

      auto cIt    = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? 1 : 0;

    }

  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

  // Preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,
                          ilu0, // preconditioner
                          1e-4, // desired residual reduction factor
                          50,   // maximum number of iterations
                          2);   // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  std::shared_ptr<VTKBasisGridFunction<FEBasis,VectorType> > xFunction
    = std::make_shared<VTKBasisGridFunction<FEBasis,VectorType> >(feBasis, x, "solution");

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(xFunction);
  vtkWriter.write("poisson-pq2");

 }
// Error handling
 catch (Exception e) {
    std::cout << e << std::endl;
 }
