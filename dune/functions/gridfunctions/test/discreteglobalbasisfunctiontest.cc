// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;


template<class B, class C>
bool checkInterpolationConsistency(B&& basis, C&& x)
{
  using Coeff = typename std::decay<C>::type;
  using Range = typename Coeff::value_type;

  // generate a discrete function
  auto f = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>
                            (basis, Dune::TypeTree::hybridTreePath(), x);

  Coeff y;
  interpolate(basis, y, f);
  for (typename std::decay_t<C>::size_type i=0; i<x.size(); ++i)
  {
    auto diff = x[i];
    diff -= y[i];
    if (diff.infinity_norm()> 1e-10)
    {
      std::cout << "Interpolation of DiscreteGlobalBasisFunction differs from original coefficient vector" << std::endl;
      return false;
    }
  }
  return true;
}


int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
  bool passed = true;

  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  // Test whether PQ1FunctionSpaceBasis.hh can be instantiated on the leaf view
  typedef GridType::LeafGridView GridView;
//  typedef PQ1NodalBasis<GridView> Basis;
  typedef PQkNodalBasis<GridView,2> Basis;

  const GridView& gridView = grid.leafGridView();
  Basis feBasis(gridView);

  using Domain = GridView::template Codim<0>::Geometry::GlobalCoordinate;

  {
    using Range = FieldVector<double,5>;
    auto f = [](const Domain& x){
      Range y;
      for (typename Range::size_type i=0; i<y.size(); ++i)
        y[i] = x[0]+i;
      return y;
    };

    std::vector<Range> x;
    interpolate(feBasis, x, f);
    passed = passed and checkInterpolationConsistency(feBasis, x);
  }


  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result. Notice that resizing is done by the interpolate method.
  std::vector<FieldVector<double,1> > x;
  auto fAnalytic = [](const Domain& x){ return x[0];};
  interpolate(feBasis, x, fAnalytic);

  passed = passed and checkInterpolationConsistency(feBasis, x);

  // generate a discrete function to evaluate the integral
  auto f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
                            (feBasis, Dune::TypeTree::hybridTreePath(), x);

  double exactIntegral = 0.5;

  std::cout << "Testing with raw DiscreteGlobalBasisFunction" << std::endl;
  passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, f, exactIntegral);

  if (passed)
    std::cout << "All tests passed" << std::endl;

  return passed ? 0: 1;

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
