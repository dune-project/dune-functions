// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/differentiablefunctionfromcallables.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;

double x_component_A(const Dune::FieldVector<double, 2> & x)
{
    return x[0];
}

template<typename Coord>
double x_component_B(const Coord & x)
{
    return x[0];
}

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
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

  std::cout << "Testing manual construction of AnalyticGridViewFunction with range type double" << std::endl;
  {
    using Range = double;

    auto f = [](const Domain& x) {return x[0];};

    AnalyticGridViewFunction<Range(Domain), GridView, decltype(f)> fGVF(f, gridView);

    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, fGVF, exactIntegral);
  }

  std::cout << "Testing deduction guide of AnalyticGridViewFunction with range type double" << std::endl;
  {
    auto f = [](const Domain& x) {return x[0];};

    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, AnalyticGridViewFunction{f, gridView}, exactIntegral);
  }

  std::cout << "Testing makeAnalyticGridViewFunction with range type double" << std::endl;
  {
    auto f = [](const Domain& x) {return x[0];};

    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, makeAnalyticGridViewFunction(f, gridView), exactIntegral);
  }

  std::cout << "Testing makeGridViewFunction with range type double" << std::endl;
  {
    auto f = [](const Domain& x) {return x[0];};

    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, makeGridViewFunction(f, gridView), exactIntegral);
  }

  std::cout << "Testing makeAnalyticGridViewFunction with free function" << std::endl;
  {
    auto f = x_component_A;
    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, makeAnalyticGridViewFunction(f, gridView), exactIntegral);
  }

  std::cout << "Testing makeAnalyticGridViewFunction with free template function" << std::endl;
  {
    auto f = x_component_B<Domain>;
    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, makeAnalyticGridViewFunction(f, gridView), exactIntegral);
  }

  std::cout << "Testing makeAnalyticGridViewFunction with derivate" << std::endl;
  {
    using Range = double;
    auto f = [](Domain){ return 1.0; };
    auto df = [](Domain){ return FieldVector<Range,2>({0.0, 0.0}); };
    auto af = makeDifferentiableFunctionFromCallables(Dune::Functions::SignatureTag<Range(Domain)>(), f, df);
    auto gf = makeAnalyticGridViewFunction(af, gridView);
    passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, gf, 1.0);

    auto ep = gridView.template begin<0>();
    {
      std::cout << "Checking evaluation of function and derivative" << std::endl;
      auto _lf = localFunction(gf);
      _lf.bind(*ep);
      auto _df = derivative(gf);
      auto _ldf = localFunction(_df);
      _ldf.bind(*ep);
      passed = passed and gf({0.0,0.0}) == 1.0;
      passed = passed and _lf({0.0,0.0}) == 1.0;
      passed = passed and _df({0.0,0.0}) == Domain(0.0);
      passed = passed and _ldf({0.0,0.0}) == Domain(0.0);
    }
    {
      std::cout << "Checking evaluation of function and derivative via Interface class" << std::endl;
      GridViewFunction<Range(Domain), GridView> _gf = gf;
      auto _lf = localFunction(_gf);
      _lf.bind(*ep);
      auto _df = derivative(_gf);
      auto _ldf = localFunction(_df);
      _ldf.bind(*ep);
      passed = passed and gf({0.0,0.0}) == 1.0;
      passed = passed and _lf({0.0,0.0}) == 1.0;
      passed = passed and _df({0.0,0.0}) == Domain(0.0);
      passed = passed and _ldf({0.0,0.0}) == Domain(0.0);
    }

    {
      std::cout << "Checking evaluation of local-function and local-derivative" << std::endl;
      auto lf = localFunction(gf);
      lf.bind(*ep);
      auto dlf = derivative(lf);
      passed = passed and gf({0.0,0.0}) == 1.0;
      passed = passed and lf({0.0,0.0}) == 1.0;
      passed = passed and dlf({0.0,0.0}) == Domain(0.0);
    }
  }

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
