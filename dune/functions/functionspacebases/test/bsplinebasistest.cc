// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/bsplinebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;



template <int dim>
void testForDimension(TestSuite& test)
{
  std::cout << "   +++++++++++  Testing on " << dim << "d grid  ++++++++++++" << std::endl;

  // Generate grid for testing
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l;
  std::fill(l.begin(), l.end(), 1.0);
  std::array<int,dim> elements;
  std::fill(elements.begin(), elements.end(), 2);
  GridType grid(l,elements);

  // Test whether function space basis can be instantiated on the leaf view
  typedef typename GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  grid.globalRefine(2);

  // Testing B-spline basis with open knot vectors
  std::vector<double> knotVector(elements[0]*4+1);
  for (size_t i=0; i<knotVector.size(); i++)
    knotVector[i] = i*l[0] / elements[0];

  // Test open knot vectors
  std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;
  for (unsigned int order : {0, 1, 2})
  {
    {
      // Check basis created via its constructor
      Functions::BSplineBasis<GridView> basis(gridView, knotVector, order);
      test.subTest(checkBasis(basis));
    }

    {
      // Check basis created via makeBasis
      using namespace Functions::BasisFactory;
      auto basis = makeBasis(gridView, bSpline(knotVector, order));
      test.subTest(checkBasis(basis));
    }

    {
      // Check whether a B-Spline basis can be combined with other bases.
      using namespace Functions::BasisFactory;
      auto basis = makeBasis(gridView,
                             power<2>(
                               bSpline(knotVector, order)
                             ));
      test.subTest(checkBasis(basis));
    }
  }

  // Testing B-spline basis with non-open knot vectors
  std::cout << "  Testing B-spline basis with non-open knot vectors" << std::endl;
  for (unsigned int order : {0, 1, 2})
  {
    {
      // Check basis created via its constructor
      Functions::BSplineBasis<GridView> bSplineBasis(gridView, knotVector, order, false);
      test.subTest(checkBasis(bSplineBasis));
    }

    {
      // Check basis created via makeBasis
      using namespace Functions::BasisFactory;
      auto basis = makeBasis(gridView, bSpline(knotVector, order, false));
      test.subTest(checkBasis(basis));
    }
  }
}


int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  testForDimension<1>(test);
  testForDimension<2>(test);
  testForDimension<3>(test);

  return test.exit();
}
