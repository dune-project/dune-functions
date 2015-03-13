// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_TEST_GRIDFUNCTIONTEST_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_TEST_GRIDFUNCTIONTEST_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>


namespace Dune {
namespace Functions {
namespace Test {



template<class GridView, class F>
double integrateGridViewFunction(const GridView& gridView, const F& f)
{
  static const int dim = GridView::dimension;

  double integral = 0;

  auto fLocal = localFunction(f);

  // Loop over elements and integrate over the function
//  for (auto it = gridView.begin<0>(); it != gridView.end<0>(); ++it)
  for (const auto& e : elements(gridView))
  {
//    const auto& e = *it;
    auto geometry = e.geometry();

    fLocal.bind(e);

    // A quadrature rule
    const auto& quad = QuadratureRules<double, dim>::rule(e.type(), 1);

    // Loop over all quadrature points
    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      auto quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      auto integrationElement = geometry.integrationElement(quadPos);

      integral += fLocal(quadPos) * quad[pt].weight() * integrationElement;
    }
    fLocal.unbind();
  }
  return integral;
}


template<class GridView, class F>
bool checkGridViewFunction(const GridView& gridView, const F& f, double exactIntegral)
{
  bool passed = true;
  double integral;

  std::cout << "Checking integration of raw function f on grid view" << std::endl;
  integral = integrateGridViewFunction(gridView, f);
  if (std::abs(integral-0.5)> 1e-10)
  {
    std::cout << "ERROR: Integral is " << integral << " but should be " << exactIntegral << std::endl;
    passed = false;
  }

  using EntitySet = typename F::EntitySet;

  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = typename std::result_of<F(Domain)>::type;

  std::cout << "Checking integration of GridFunction<Range(Domain), EntitySet>(f) on grid view" << std::endl;
  GridFunction<Range(Domain), EntitySet> f2 = f;
  integral = integrateGridViewFunction(gridView, f2);
  if (std::abs(integral-exactIntegral)> 1e-10)
  {
    std::cout << "ERROR: Integral is " << integral << " but should be " << exactIntegral << std::endl;
    passed = false;
  }

  std::cout << "Checking integration of GridViewFunction<Range(Domain), GridView>(f) on grid view" << std::endl;
  GridViewFunction<Range(Domain), GridView> f3 = f;
  integral = integrateGridViewFunction(gridView, f3);
  if (std::abs(integral-0.5)> 1e-10)
  {
    std::cout << "ERROR: Integral is " << integral << " but should be " << exactIntegral << std::endl;
    passed = false;
  }

  std::cout << "Checking integration of makeGridFunction(f) on grid view" << std::endl;
  auto f4 = makeGridViewFunction(f, gridView);
  integral = integrateGridViewFunction(gridView, f4);
  if (std::abs(integral-0.5)> 1e-10)
  {
    std::cout << "ERROR: Integral is " << integral << " but should be " << exactIntegral << std::endl;
    passed = false;
  }

  std::cout << "Checking integration of GridViewFunction<Range(Domain), GridView>(makeGridFunction(f)) on grid view" << std::endl;
  GridViewFunction<Range(Domain), GridView> f5 = f4;
  integral = integrateGridViewFunction(gridView, f5);
  if (std::abs(integral-0.5)> 1e-10)
  {
    std::cout << "ERROR: Integral is " << integral << " but should be " << exactIntegral << std::endl;
    passed = false;
  }

  return passed;
}




} // namespace Test
} // namespace Functions
} // namespace Dune




#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_TEST_GRIDFUNCTIONTEST_HH
