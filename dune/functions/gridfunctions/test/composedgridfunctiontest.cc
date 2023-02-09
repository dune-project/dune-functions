// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/istl/bvector.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/composedgridfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  TestSuite suite;

  // Generate grid for testing
  const int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  auto l = FieldVector<double,dim>{1, 1};
  auto elements = std::array<int,dim>{{10, 10}};
  auto grid = Grid(l, elements);
  auto gridView = grid.leafGridView();

  using namespace Functions::BasisBuilder;

  {
    using Range = FieldVector<double,2>;

    // Inner test function f is a polynomial of degree 2.
    auto f = [](const auto& x){
      Range y;
      for (typename Range::size_type i=0; i<y.size(); ++i)
        y[i] = (x[i]+i)*x[i];
      return y;
    };

    // Outer test function f is a polynomial of degree 2.
    auto g = [](auto v0, auto v1) {
      return (v0*v0) + (v1*v1);
    };

    // Coposition is a polynomial of degree 4.
    auto gf = [&](auto x) {
      return g(f(x)[0], f(x)[1]);
    };

    // Compute integral (order 4 is sufficient).
    auto integral = integrateGridViewFunction(gridView, makeAnalyticGridViewFunction(gf, gridView), 4);

    // Create 2nd order polynomial basis
    auto basis = makeBasis(gridView,power<2>(lagrange<2>()));
    Dune::BlockVector<Range> c;

    // Interpolate f wrt basis.
    interpolate(basis, c, f);

    // Create grid functions from coefficients of individual
    // components. Each is the piecewise Q2 interpolation of f_i.
    // They should be exact, since f is P2.
    auto basis_0 = Dune::Functions::subspaceBasis(basis, 0);
    auto basis_1 = Dune::Functions::subspaceBasis(basis, 1);
    auto f0_gridfunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis_0, c);
    auto f1_gridfunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis_1, c);

    // Now compose on top of individual component GridFunctions and check.
    {
      // Check with capture by std::ref
      auto gf_gridfunction = makeComposedGridFunction(g, std::ref(f0_gridfunction), std::ref(f1_gridfunction));
      auto gf_gridfunction_= ComposedGridFunction(g, std::ref(f0_gridfunction), std::ref(f1_gridfunction));
      static_assert(std::is_same_v<decltype(gf_gridfunction), decltype(gf_gridfunction_)>);
      suite.check(
          checkGridViewFunction(gridView, gf_gridfunction, integral, 4),
          "Check if ComposedGridFunction has correct integral (capture with std::ref)");
    }
    {
      // Check with capture by std::cref
      auto gf_gridfunction = makeComposedGridFunction(g, std::cref(f0_gridfunction), std::cref(f1_gridfunction));
      auto gf_gridfunction_= ComposedGridFunction(g, std::cref(f0_gridfunction), std::cref(f1_gridfunction));
      static_assert(std::is_same_v<decltype(gf_gridfunction), decltype(gf_gridfunction_)>);
      suite.check(
          checkGridViewFunction(gridView, gf_gridfunction, integral, 4),
          "Check if ComposedGridFunction has correct integral (capture with std::cref)");
    }
    {
      // Check with capture by value
      auto gf_gridfunction = makeComposedGridFunction(g, f0_gridfunction, f1_gridfunction);
      auto gf_gridfunction_= ComposedGridFunction(g, f0_gridfunction, f1_gridfunction);
      static_assert(std::is_same_v<decltype(gf_gridfunction), decltype(gf_gridfunction_)>);
      suite.check(
          checkGridViewFunction(gridView, gf_gridfunction, integral, 4),
          "Check if ComposedGridFunction has correct integral (capture by value)");
    }
  }


  return suite.exit();

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
