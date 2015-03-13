// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/functionspacebases/pq2nodalbasis.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>


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

  // Test whether PQ1FunctionSpaceBasis.hh can be instantiated on the leaf view
  typedef GridType::LeafGridView GridView;
  typedef PQ1NodalBasis<GridView> Basis;
//  typedef PQ2NodalBasis<GridView> Basis;

  const GridView& gridView = grid.leafGridView();
  Basis feBasis(gridView);
  auto indexSet = feBasis.indexSet();

  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result
  std::vector<double> x(indexSet.size());

  // TODO: Implement interpolation properly using the global basis.
  for (auto it = gridView.begin<dim>(); it != gridView.end<dim>(); ++it)
    x[gridView.indexSet().index(*it)] = it->geometry().corner(0)[0];

  // generate a discrete function to evaluate the integral
  Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasis),decltype(x)> f(feBasis,x);



  double exactIntegral = 0.5;
  bool passed = true;

  passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, f, exactIntegral);

  using F = decltype(f);

  using Range = F::Range;
  using Domain = F::Domain;

  auto fAnalytic = [](const Domain& x) {return x[0];};
  AnalyticGridViewFunction<Range(Domain), GridView, decltype(fAnalytic)> f2(fAnalytic, gridView);

  passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, f2, exactIntegral);






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
