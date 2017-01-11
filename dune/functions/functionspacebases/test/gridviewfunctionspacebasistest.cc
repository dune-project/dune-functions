// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/bsplinebasis.hh>

using namespace Dune;
using namespace Dune::Functions;

/** \param disableInterpolate Not all bases can correctly represent the linear function
 *      that we use to test whether integrating over a given function gives the correct
 *      value (e.g. basis without dofs on the boundary).  Therefore we allow to disable this test.
 *  \param disabledLocalTests Allows to disable certain local tests (see dune-localfunctions/dune/localfunctions/test/test-localfe.hh)
 */
template <typename Basis>
void testScalarBasisConst(const Basis& feBasis,
                     bool isPartitionOfUnity,
                     bool disableInterpolate = false,
                     char disabledLocalTests = DisableNone)
{
  static const int dim = Basis::GridView::dimension;

  //////////////////////////////////////////////////////////////////////////////////////
  //  Run the dune-localfunctions test for the LocalFiniteElement of each grid element
  //////////////////////////////////////////////////////////////////////////////////////

  typedef typename Basis::GridView GridView;
  GridView gridView = feBasis.gridView();

  typename Basis::LocalView localView(feBasis);


  // Test the LocalFiniteElement
  for (auto it = gridView.template begin<0>(); it!=gridView.template end<0>(); ++it)
  {
    // Bind the local FE basis view to the current element
    localView.bind(*it);

    // The general LocalFiniteElement unit test from dune/localfunctions/test/test-localfe.hh
    const auto& lFE = localView.tree().finiteElement();
//    testFE(lFE, disabledLocalTests);

    if (lFE.size() != localView.size())
      DUNE_THROW(Exception, "Size of leaf node and finite element do not coincide");
  }

  /////////////////////////////////////////////////////////////////////////
  //  Make sure the basis is a partition of unity
  /////////////////////////////////////////////////////////////////////////
  if (isPartitionOfUnity)
  {
    for(const auto& e : elements(gridView))
    {
      // Bind the local FE basis view to the current element
      localView.bind(e);

      const auto& lFE = localView.tree().finiteElement();

      const QuadratureRule<double,dim>& quad = QuadratureRules<double,dim>::rule(e.type(), 3);
      std::vector<FieldVector<double,1> > values;
      for (size_t i=0; i<quad.size(); i++)
      {
        lFE.localBasis().evaluateFunction(quad[i].position(), values);
        double sum = std::accumulate(values.begin(), values.end(), 0.0);

        if (std::abs(sum-1.0) > 1e-5)
          DUNE_THROW(Exception, "Basis is no partition of unity, even though it is supposed to be! Error occurred for geometry type: " << e.type());
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  //  Make sure the basis does not give out constant zero shape functions
  ////////////////////////////////////////////////////////////////////////////
  for(const auto& e : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(e);

    const auto& lFE = localView.tree().finiteElement();

    const QuadratureRule<double,dim>& quad = QuadratureRules<double,dim>::rule(e.type(), 3);
    std::vector<FieldVector<double,1> > values;
    std::vector<double> sumOfAbsValues(lFE.size());
    std::fill(sumOfAbsValues.begin(), sumOfAbsValues.end(), 0.0);
    for (size_t i=0; i<quad.size(); i++)
    {
      lFE.localBasis().evaluateFunction(quad[i].position(), values);
      // Sum the absolute values for each quadrature point
      for (size_t j=0; j<values.size(); j++)
        sumOfAbsValues[j] += std::fabs(values[j][0]);
    }

    for (size_t i=0; i<values.size(); i++)
      if ( sumOfAbsValues[i] <= 1e-5)
        DUNE_THROW(Exception, "Basis gives out a constant-zero shape function!");
  }




  // Check whether the basis exports a type 'MultiIndex'
  typedef typename Basis::MultiIndex MultiIndex;

  // And this type must be indexable
  static_assert(is_indexable<MultiIndex>(), "MultiIndex must support operator[]");

  ///////////////////////////////////////////////////////////////////////////////////
  //  Check whether the global indices are in the correct range,
  //  and whether each global index appears at least once.
  ///////////////////////////////////////////////////////////////////////////////////

  std::vector<bool> seen(feBasis.size());
  std::fill(seen.begin(), seen.end(), false);

  auto localIndexSet = feBasis.localIndexSet();

  // Loop over all leaf elements
  for (auto it = gridView.template begin<0>(); it!=gridView.template end<0>(); ++it)
  {
    // Bind the local FE basis view to the current element
    localView.bind(*it);
    localIndexSet.bind(localView);

    for (size_t i=0; i<localView.tree().size(); i++)
    {
      if (localIndexSet.index(i)[0] < 0)
        DUNE_THROW(Exception, "Index is negative, which is not allowed");

      if (localIndexSet.index(i)[0] >= seen.size())
        DUNE_THROW(Exception, "Local index " << i
                           << " is mapped to global index " << localIndexSet.index(i)
                           << ", which is larger than allowed");

      seen[localIndexSet.index(i)[0]] = true;
    }
  }

  for (size_t i=0; i<seen.size(); i++)
    if (! seen[i])
      DUNE_THROW(Exception, "Index [" << i << "] does not exist as global basis vector");

  //////////////////////////////////////////////////////////////////////////////////////////
  // Interpolate the function f(x,y) = x wrt the basis, and check whether we get
  // the expected integral.
  //////////////////////////////////////////////////////////////////////////////////////////

  std::vector<double> x(feBasis.size());
  if (! disableInterpolate)
    interpolate(feBasis, x, [](FieldVector<double,dim> x){ return x[0]; });
  else  // dummy values
    std::fill(x.begin(), x.end(), 0.5);

  // Objects required in the local context
  auto localIndexSet2 = feBasis.localIndexSet();
  std::vector<double> localCoefficients(localView.maxSize());

  // Loop over elements and integrate over the function
  double integral = 0;
  for (const auto& element : elements(gridView))
  {
    localView.bind(element);
    localIndexSet.bind(localView);
    localIndexSet2.bind(localView);

    // paranoia checks
    assert(localView.size() == localIndexSet.size());
    assert(&(localView.globalBasis()) == &(feBasis));
    assert(&(localIndexSet.localView()) == &(localView));

    assert(localIndexSet.size() == localIndexSet2.size());
    for (size_t i=0; i<localIndexSet.size(); i++)
      assert(localIndexSet.index(i) == localIndexSet2.index(i));

    // copy data from global vector
    localCoefficients.resize(localIndexSet.size());
    for (size_t i=0; i<localIndexSet.size(); i++)
      localCoefficients[i] = x[localIndexSet.index(i)[0]];

    // get access to the finite element
    typedef typename Basis::LocalView::Tree Tree;
    const Tree& tree = localView.tree();

    auto& localFiniteElement = tree.finiteElement();

    // we have a flat tree...
    assert(localView.size() == tree.size());
    assert(localView.size() == tree.finiteElement().localBasis().size());

    // A quadrature rule
    // BUG: I need more than just the order given by the localBasis to make the integral come out precisely
    const auto& quad = QuadratureRules<double, dim>::rule(element.type(), 1+tree.finiteElement().localBasis().order());

    // Loop over all quadrature points
    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = element.geometry().integrationElement(quadPos);

      // Evaluate all shape function values at this point
      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

      // Actually compute the vector entries
      for (size_t i=0; i<localFiniteElement.localBasis().size(); i++)
      {
        integral += localCoefficients[tree.localIndex(i)] * shapeFunctionValues[i] * quad[pt].weight() * integrationElement;
      }
    }

    // unbind
    localIndexSet.unbind();
    localView.unbind();
  }

  if (!disableInterpolate && std::abs(integral-0.5) > 1e-10)
    DUNE_THROW(Dune::Exception, "Error: integral value is wrong!");
}

template <typename Basis, typename GV>
void testScalarBasis(Basis& feBasis, GV gv,
                     bool isPartitionOfUnity,
                     bool disableInterpolate = false,
                     char disabledLocalTests = DisableNone)
{
  feBasis.update(gv);
  testScalarBasisConst(feBasis, isPartitionOfUnity, disableInterpolate, disabledLocalTests);
}

template <int dim>
void testOnStructuredGrid()
{
  std::cout << "   +++++++++++  Testing on structured " << dim << "d grids  ++++++++++++" << std::endl;

  // Generate grid for testing
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l;
  std::fill(l.begin(), l.end(), 1.0);
  std::array<int,dim> elements;
  std::fill(elements.begin(), elements.end(), 2);
  GridType grid(l,elements);

  // Test whether PQ1FunctionSpaceBasis.hh can be instantiated on the leaf view
  typedef typename GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  PQ1NodalBasis<GridView> pq1Basis(gridView);
  PQkNodalBasis<GridView, 3> pq3Basis(gridView);
  PQkNodalBasis<GridView, 4> pq4Basis(gridView);
  PQkNodalBasis<GridView, 0> pq0Basis(gridView);
  LagrangeDGBasis<GridView, 1> lagrangeDG1Basis(gridView);
  LagrangeDGBasis<GridView, 2> lagrangeDG2Basis(gridView);
  LagrangeDGBasis<GridView, 3> lagrangeDG3Basis(gridView);

  grid.globalRefine(2);

  // Test PQ1NodalBasis
  testScalarBasis(pq1Basis, gridView, true);

  // Test PQkNodalBasis for k==3
  if (dim<3) // Currently not implemented for dim >= 3
    testScalarBasis(pq3Basis, gridView, true);

  // Test PQkNodalBasis for k==4
  if (dim<3) // Currently not implemented for dim >= 3
    testScalarBasis(pq4Basis, gridView, true);

  // Test PQkNodalBasis for the corner case k == 0
  testScalarBasis(pq0Basis, gridView, true);

  // Test LagrangeDGBasis for k==1
  testScalarBasis(lagrangeDG1Basis, gridView, true);

  // Test LagrangeDGBasis for k==2
  testScalarBasis(lagrangeDG2Basis, gridView, true);

  // Test LagrangeDGBasis for k==3
  testScalarBasis(lagrangeDG3Basis, gridView, true);

  // Testing B-spline basis with open knot vectors
  std::vector<double> knotVector(elements[0]*4+1);
  for (size_t i=0; i<knotVector.size(); i++)
    knotVector[i] = i*l[0] / elements[0];

  // Test open knot vectors
  std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;
  for (unsigned int order : {0, 1, 2})
  {
    std::cout << "   order: " << order << std::endl;

    BSplineBasis<GridView> bSplineBasis(gridView, knotVector, order);
    testScalarBasisConst(bSplineBasis,
                    true,
                    true,      // Don't interpolate a given function and try to integrate over it
                    DisableLocalInterpolation);
  }

  // Testing B-spline basis with non-open knot vectors
  std::cout << "  Testing B-spline basis with non-open knot vectors" << std::endl;
  for (unsigned int order : {0, 1, 2})
  {
    std::cout << "   order: " << order << std::endl;

    BSplineBasis<GridView> bSplineBasis(gridView, knotVector, order, false);
    testScalarBasisConst(bSplineBasis,
                    order==0,   // only zero-order B-splines for a partition of unity
                    true,      // Don't interpolate a given function and try to integrate over it
                    DisableLocalInterpolation);
  }
}

template <int dim>
void testOnHybridGrid()
{
#if ! HAVE_UG
  std::cout << "   ---- Skipping test on hybrid " << dim << "d grid -- UGGrid is not installed ----" << std::endl;
#else
  std::cout << "   +++++++++++  Testing on hybrid " << dim << "d grid  ++++++++++++" << std::endl;

  // Generate grid for testing
  typedef UGGrid<dim> GridType;

  const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

  std::string filename = path + "hybrid-testgrid-" + std::to_string(dim) + "d.msh";

  std::shared_ptr<GridType> grid(GmshReader<GridType>::read(filename));

  // Test whether function space basis can be instantiated on the leaf view
  typedef typename GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  // Disable the global interpolation test for 3d grids.
  // Global interpolation should work for all grids and spaces, but the hard-coded
  // reference value is wrong.
  bool disableInterpolate = (dim==3);

  // Test PQ1NodalBasis -- dedicated implementation
  PQ1NodalBasis<GridView> pq1DedicatedBasis(gridView);
  testScalarBasisConst(pq1DedicatedBasis, true, disableInterpolate);

  // Test PQ1NodalBasis -- generic basis
  PQkNodalBasis<GridView, 1> pq1Basis(gridView);
  testScalarBasisConst(pq1Basis, true, disableInterpolate);

  // Test PQkNodalBasis for k==3
  PQkNodalBasis<GridView, 3> pq3Basis(gridView);
  if (dim<3) // Currently not implemented for dim >= 3
    testScalarBasisConst(pq3Basis, true, disableInterpolate);

  // Test PQkNodalBasis for k==4
  PQkNodalBasis<GridView, 4> pq4Basis(gridView);
  if (dim<3) // Currently not implemented for dim >= 3
    testScalarBasisConst(pq4Basis, true, disableInterpolate);

  // Test LagrangeDGBasis for k==1
  LagrangeDGBasis<GridView, 1> lagrangeDG1Basis(gridView);
  testScalarBasisConst(lagrangeDG1Basis, true, disableInterpolate);

  // Test LagrangeDGBasis for k==2
  // \todo Enable these tests once pyramid element of order two is bug free
//  LagrangeDGBasis<GridView, 2> lagrangeDG2Basis(gridView);
//  testScalarBasisConst(lagrangeDG2Basis, true, disableInterpolate);

  // Test LagrangeDGBasis for k==3
  // \todo Enable these tests once pyramid element of order three is implemented
//  LagrangeDGBasis<GridView, 3> lagrangeDG3Basis(gridView);
//  testScalarBasisConst(lagrangeDG3Basis, true, disableInterpolate);
#endif
}



int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  testOnStructuredGrid<1>();
  testOnStructuredGrid<2>();
  testOnStructuredGrid<3>();

  testOnHybridGrid<2>();
  testOnHybridGrid<3>();

  return 0;

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
