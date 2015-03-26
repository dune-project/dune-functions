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


#define CHECK(B,M) ( [&]() { if (not(B)) std::cout << "TEST FAILURE:" << M << std::endl; return B;}() )

template<class GridView, class F>
double integrateGridViewFunction(const GridView& gridView, const F& f, unsigned int quadOrder)
{
  static const int dim = GridView::dimension;

  double integral = 0;

  auto fLocal = localFunction(f);

  // Loop over elements and integrate over the function
  for (const auto& e : elements(gridView))
  {
    auto geometry = e.geometry();

    fLocal.bind(e);

    // A quadrature rule
    const auto& quad = QuadratureRules<double, dim>::rule(e.type(), quadOrder);

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
bool checkGridViewFunction(const GridView& gridView, const F& f, double exactIntegral, unsigned int quadOrder=1)
{
  bool passed = true;
  double integral;

  std::cout << "Checking integration of raw function f on grid view" << std::endl;
  integral = integrateGridViewFunction(gridView, f, quadOrder);
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
  integral = integrateGridViewFunction(gridView, f2, quadOrder);
  passed = CHECK(std::abs(integral-exactIntegral) < 1e-10, "Integral is " << integral << " but should be " << exactIntegral);

  std::cout << "Checking integration of GridViewFunction<Range(Domain), GridView>(f) on grid view" << std::endl;
  GridViewFunction<Range(Domain), GridView> f3 = f;
  integral = integrateGridViewFunction(gridView, f3, quadOrder);
  passed = CHECK(std::abs(integral-exactIntegral) < 1e-10, "Integral is " << integral << " but should be " << exactIntegral);

  std::cout << "Checking integration of makeGridFunction(f) on grid view" << std::endl;
  auto f4 = makeGridViewFunction(f, gridView);
  integral = integrateGridViewFunction(gridView, f4, quadOrder);
  passed = CHECK(std::abs(integral-exactIntegral) < 1e-10, "Integral is " << integral << " but should be " << exactIntegral);

  std::cout << "Checking integration of GridViewFunction<Range(Domain), GridView>(makeGridFunction(f)) on grid view" << std::endl;
  GridViewFunction<Range(Domain), GridView> f5 = f4;
  integral = integrateGridViewFunction(gridView, f5, quadOrder);
  passed = CHECK(std::abs(integral-exactIntegral) < 1e-10, "Integral is " << integral << " but should be " << exactIntegral);

  std::cout << "Checking integration of default constructed and assigned GridViewFunction<Range(Domain), GridView> on grid view" << std::endl;
  GridViewFunction<Range(Domain), GridView> f6;
  f6 = f5;
  integral = integrateGridViewFunction(gridView, f6, quadOrder);
  passed = CHECK(std::abs(integral-exactIntegral) < 1e-10, "Integral is " << integral << " but should be " << exactIntegral);

  std::cout << "Checking integration of reassigned GridViewFunction<Range(Domain), GridView> on grid view" << std::endl;
  f6 = f3;
  integral = integrateGridViewFunction(gridView, f6, quadOrder);
  passed = CHECK(std::abs(integral-exactIntegral) < 1e-10, "Integral is " << integral << " but should be " << exactIntegral);

  return passed;
}




} // namespace Test
} // namespace Functions
} // namespace Dune




#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_TEST_GRIDFUNCTIONTEST_HH
