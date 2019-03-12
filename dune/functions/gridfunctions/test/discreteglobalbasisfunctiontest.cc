// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;


int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  // Test whether LagrangeBasis can be instantiated on the leaf view
  typedef GridType::LeafGridView GridView;
  typedef LagrangeBasis<GridView,2> Basis;

  const GridView& gridView = grid.leafGridView();
  Basis feBasis(gridView);

  using Domain = GridView::template Codim<0>::Geometry::GlobalCoordinate;

  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result. Notice that resizing is done by the interpolate method.
  std::vector<FieldVector<double,1> > x(feBasis.size());
  std::fill(x.begin(), x.end(), 0);

  // generate a discrete function
  auto f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
                            (feBasis, Dune::TypeTree::hybridTreePath(), x);

  using F = decltype(f);
  using EntitySet = typename F::EntitySet;
  using Range = typename std::result_of<F(Domain)>::type;
  GridViewFunction<Range(Domain), GridView> f3 = f;
  checkTrue(Functions::Concept::isGridFunction<decltype(f3), Range(Domain), EntitySet>(), "Function does not model GridFunction concept");

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
