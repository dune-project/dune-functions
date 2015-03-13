// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;

int main (int argc, char* argv[]) try
{
  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  using GridView = typename GridType::LeafGridView;

  const GridView& gridView = grid.leafGridView();


  double exactIntegral = 0.5;
  bool passed = true;

  using Domain = GridView::Codim<0>::Geometry::GlobalCoordinate;

  std::cout << "Testing with range type double" << std::endl;
  {
    using Range = double;

    auto f = [](const Domain& x) {return x[0];};

    AnalyticGridViewFunction<Range(Domain), GridView, decltype(f)> fGVF(f, gridView);

    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, fGVF, exactIntegral);
  }






  if (passed)
    std::cout << "All tests passed" << std::endl;

  return passed ? 0: 1;

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
}
