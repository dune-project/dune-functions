// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <cmath>
#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/common/functionconcepts.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/coarsefunctiononfinegridview.hh>
#include <dune/functions/gridfunctions/finefunctiononcoarsegridview.hh>



template<class X>
requires Dune::IsNumber<X>::value
auto rangeDist(X x1, X x2)
{
  return std::abs(x1-x2);
}

template<class F, int n>
auto rangeDist(Dune::FieldVector<F, n> x1, Dune::FieldVector<F, n> x2)
{
  x1 -= x2;
  return x1.infinity_norm();
}

template<class F, int n, int m>
auto rangeDist(Dune::FieldMatrix<F, n, m> x1, Dune::FieldMatrix<F, n, m> x2)
{
  x1 -= x2;
  return x1.infinity_norm();
}



// This is used as a drop-in replacement for Dune::QuadratureRules
// except for pyramids. For the latter, the standard Dune rules provide
// points on the diagonal. Unfortunately the compile time order Lagrange
// elements on pyramids are piecewise polynomial with gradients jumping
// across the diagonal. This utility defines pyramid rules by splitting
// them into two tetrahedra. Since the pyramid order-k shape functions
// are bi-order-k polynomials double the order for those.
template<class ctype, int dim>
struct QuadratureRules
{

  using QuadratureRule = Dune::QuadratureRule<ctype, dim>;
  using QuadraturePoint = Dune::QuadraturePoint<ctype, dim>;
  using Coordinate = typename QuadraturePoint::Vector;

  static const auto& rule(const Dune::GeometryType& type, size_t order)
  {
    static auto pyramidRules = std::vector<QuadratureRule>();
    if (type == Dune::GeometryTypes::pyramid)
    {
      if (pyramidRules.size() < order)
        for(auto k : Dune::range(pyramidRules.size(), std::size_t(order+1)))
        {
          const auto& simplexRule = Dune::QuadratureRules<double, dim>::rule(Dune::GeometryTypes::tetrahedron, 2*k);
          auto pyramidRule = QuadratureRule();
          for(auto [position, weight] : simplexRule)
          {
            pyramidRule.emplace_back(Coordinate({position[0], position[0]+position[1], position[2]}), weight);
            pyramidRule.emplace_back(Coordinate({position[0]+position[1], position[1], position[2]}), weight);
          }
          pyramidRules.push_back(pyramidRule);
        }
      return pyramidRules[order];
    }
    else
      return Dune::QuadratureRules<double, dim>::rule(type, order);
  }
};



// Check difference of local functions on element
template<class Element, class F, class G>
double maxNormDistanceOnElement(const Element& element, const F& f_local, const G& g_local, bool checkCorners, std::size_t sampleQuadratureOrder = 5)
{
  constexpr auto dimension = Element::Geometry::mydimension;
  const auto& re = referenceElement(element);
  double diff = 0.0;
  const auto& quad = QuadratureRules<double, dimension>::rule(element.type(), sampleQuadratureOrder);
  for (const auto& [x, weight] : quad)
  {
    auto f_x = f_local(x);
    auto g_x = g_local(x);
    diff = std::max(diff, rangeDist(f_x, g_x));
  }
  if (checkCorners)
    for (const auto& k : Dune::range(re.size(dimension)))
    {
      auto x = re.position(k, dimension);
      auto f_x = f_local(x);
      auto g_x = g_local(x);
      diff = std::max(diff, rangeDist(f_x, g_x));
    }
  return diff;
}

// Check difference of grid functions on grid view
template<class GridView, class F, class G>
double maxNormDistanceOnGridView(const GridView& gridView, const F& f, const G& g, bool checkVertices, std::size_t sampleQuadratureOrder = 5)
{
  double diff = 0.0;

  auto f_local = localFunction(f);
  auto g_local = localFunction(g);
  for(const auto& element : elements(gridView))
  {
    f_local.bind(element);
    g_local.bind(element);
    diff = std::max(diff, maxNormDistanceOnElement(element, f_local, g_local, checkVertices, sampleQuadratureOrder));
  }

  return diff;
}

template<class GridView, class F, class G>
double maxNormDerivativeDistanceOnGridView(const GridView& gridView, const F& f, const G& g, bool checkVertices, std::size_t sampleQuadratureOrder = 5)
{
  double diff = 0.0;

  // First check by obtaining the derivative globally
  diff = maxNormDistanceOnGridView(gridView, derivative(f), derivative(g), checkVertices, sampleQuadratureOrder);

  // Now check by obtaining the derivative locally
  auto f_local = localFunction(f);
  auto g_local = localFunction(g);
  for(const auto& element : elements(gridView))
  {
    f_local.bind(element);
    g_local.bind(element);
    auto df_local = derivative(f_local);
    auto dg_local = derivative(g_local);
    diff = std::max(diff, maxNormDistanceOnElement(element, df_local, dg_local, checkVertices, sampleQuadratureOrder));
  }

  return diff;
}


template<class Grid>
Dune::TestSuite checkOnGrid(const Grid& grid, std::string name)
{
  auto testSuite = Dune::TestSuite(name);

  using Dune::Functions::Concept::isGridViewFunction;
  using Element = typename Grid::template Codim<0>::Entity;
  using Domain = typename Element::Geometry::GlobalCoordinate;
  using Range = Dune::FieldVector<double,2>;

  double TOL = 1e-14;
  double Gradient_TOL = TOL*std::pow(2., double(grid.maxLevel()));

  std::cout << "Checking grid with " << (grid.maxLevel()+1) << " levels with tolerance " << TOL << std::endl;

  auto g = [](auto x) {
    auto y = Range();
    y[0] = 1.0;
    for(auto xi : x)
      y[0] *= std::sin(xi);
    y[1] = 1.0;
    for(auto xi : x)
      y[1] *= std::cos(xi);
    return y;
  };

  using namespace Dune::Functions::BasisFactory;
  auto preBasisFactory = power<2>(lagrange<2>(), flatInterleaved());

  auto gridView_f = grid.leafGridView();

  auto basis_f = makeBasis(gridView_f, preBasisFactory);
  std::vector<double> coeff_f;
  Dune::Functions::interpolate(basis_f, coeff_f, g);
  auto g_f = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis_f, coeff_f);

  // Check with any level grid view as coarse grid view
  for(auto level : Dune::range(grid.maxLevel()+1))
  {
    auto gridView_c = grid.levelGridView(level);

    // ********************************************************************************
    // First check different combinations of mapping the fine interpolation across grid views.
    // ********************************************************************************

    // Map fine function to coarse level grid view
    auto g_fc = Dune::Functions::FineFunctionOnCoarseGridView(g_f, gridView_c);
    testSuite.check(isGridViewFunction<decltype(g_fc), Range(Domain), decltype(gridView_c)>());

    // Remap coarse function to fine leaf grid view
    auto g_fcf = Dune::Functions::CoarseFunctionOnFineGridView(g_fc, gridView_f);
    testSuite.check(isGridViewFunction<decltype(g_fcf), Range(Domain), decltype(gridView_f)>());

    // Compare values on fine level
    {
      auto dist = maxNormDistanceOnGridView(gridView_f, g_fcf, g_f, true);
      testSuite.check(dist < TOL)
        << "Values of mapped fine->coarse->fine function differ from fine function by " << dist;
    }

    // Compare jacobians on fine level
    // Since the sample points are in the interiour of fine elements
    // they are in the interiour of coarse elements, too. Hence we
    // should not get problems due to the fact that the derivative
    // is discontinuous across fine elements. However, we should not
    // check the element corners.
    {
      auto dist = maxNormDerivativeDistanceOnGridView(gridView_f, g_fcf, g_f, false);
      testSuite.check(dist < Gradient_TOL)
        << "Jacobians of mapped fine->coarse->fine function differ from fine function by " << dist;
    }

    // Remap fine function to coarse level grid view
    auto g_fcfc = Dune::Functions::FineFunctionOnCoarseGridView(g_fcf, gridView_c);
    testSuite.check(isGridViewFunction<decltype(g_fcfc), Range(Domain), decltype(gridView_c)>());

    // Compare values on coarse level
    {
      auto dist = maxNormDistanceOnGridView(gridView_c, g_fcfc, g_fc, true);
      testSuite.check(dist < TOL)
        << "Values of mapped fine->coarse->fine->coarse function differ from fine->coarse function by " << dist;
    }

    // Do not compare jacobians on coarse level

    // We do not compare the derivatives of g_fc anfd g_fcfc on coarse
    // elements. Since the derivatives are discontinuous across fine
    // elements, the sample points within coarse elements may be on
    // a discontinuity such that both evaluate in different fine
    // elements leading to different results. Since the interpolated
    // function is C^2, the error will scale by a power of the fine
    // mesh size.
    //

    // ********************************************************************************
    // Now check different combinations of mapping the coarse interpolation across grid views.
    // ********************************************************************************

    auto basis_c = makeBasis(gridView_c, preBasisFactory);
    std::vector<double> coeff_c;
    Dune::Functions::interpolate(basis_c, coeff_c, g);
    auto g_c = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis_c, coeff_c);

    // Map coarse function to fine leaf grid view
    auto g_cf = Dune::Functions::CoarseFunctionOnFineGridView(g_c, gridView_f);
    testSuite.check(isGridViewFunction<decltype(g_cf), Range(Domain), decltype(gridView_f)>());

    // Remap fine function to coarse level grid view
    auto g_cfc = Dune::Functions::FineFunctionOnCoarseGridView(g_cf, gridView_c);
    testSuite.check(isGridViewFunction<decltype(g_cfc), Range(Domain), decltype(gridView_c)>());

    // Compare values on coarse level
    {
      auto dist = maxNormDistanceOnGridView(gridView_c, g_cfc, g_c, true);
      testSuite.check(dist < TOL)
        << "Values of mapped coarse->fine->coarse function differ from coarse function by " << dist;
    }

    // Compare jacobians on coarse level
    // Since the coarse interpolates derivative is continuous within
    // the coarse elements, this is safe. Since we only search children
    // of the coares element, we can even check the corners.
    {
      auto dist = maxNormDerivativeDistanceOnGridView(gridView_c, g_cfc, g_c, true);
      testSuite.check(dist < Gradient_TOL)
        << "Jacobians of mapped coarse->fine->coarse function differ from coarse function by " << dist;
    }

  }

  return testSuite;
}

template<class Grid>
void refineLocalNearOrigin(Grid& grid, std::size_t refinements)
{
  constexpr auto dimension = Grid::dimension;
  auto origin = Dune::FieldVector<double, dimension>();
  for([[maybe_unused]] auto k : Dune::range(refinements))
  {
    for(const auto& element : elements(grid.leafGridView()))
    {
      const auto& re = referenceElement(element);
      const auto& geometry = element.geometry();
      if (re.checkInside(geometry.local(origin)))
        grid.mark(1, element);
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();
  }
}

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite;

  {
    using Grid = Dune::YaspGrid<1>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({{0}}, {{1}}, {{2}});
    grid->globalRefine(2);
    testSuite.subTest(checkOnGrid(*grid, "YaspGrid<1>"));
  }

  {
    using Grid = Dune::YaspGrid<2>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({{0,0}}, {{1,1}}, {{2,2}});
    grid->globalRefine(2);
    testSuite.subTest(checkOnGrid(*grid, "YaspGrid<2>"));
  }

  {
    using Grid = Dune::YaspGrid<3>;
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({{0,0,0}}, {{1,1,1}}, {{2,2,2}});
    grid->globalRefine(2);
    testSuite.subTest(checkOnGrid(*grid, "YaspGrid<3>"));
  }

  // We use UGGrid with RefinementType=copy to ensure that intermediate levels cover the whole domain

  {
    using Grid = Dune::UGGrid<2>;
    auto grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid({{0,0}}, {{1,1}}, {{1,1}});
    grid->setRefinementType(Grid::RefinementType::COPY);
    refineLocalNearOrigin(*grid, 15);
    testSuite.subTest(checkOnGrid(*grid, "UGGrid<2> (triangles)"));
  }

  {
    using Grid = Dune::UGGrid<2>;
    auto factory = Dune::GridFactory<Grid>();
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({1,0});
    factory.insertVertex({2,2});
    factory.insertElement(Dune::GeometryTypes::cube(2), {0,1,2,3});
    auto grid = factory.createGrid();
    grid->setRefinementType(Grid::RefinementType::COPY);
    refineLocalNearOrigin(*grid, 15);
    testSuite.subTest(checkOnGrid(*grid, "UGGrid<2> (nonaffine rectangles)"));
  }

  {
    using Grid = Dune::UGGrid<3>;
    auto grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid({{0,0,0}}, {{1,1,1}}, {{1,1,1}});
    grid->setRefinementType(Grid::RefinementType::COPY);
    refineLocalNearOrigin(*grid, 4);
    testSuite.subTest(checkOnGrid(*grid, "UGGrid<2> (triangles)"));
  }

  {
    using Grid = Dune::UGGrid<3>;
    auto factory = Dune::GridFactory<Grid>();
    factory.insertVertex({0,0,0});
    factory.insertVertex({1,0,0});
    factory.insertVertex({0,1,0});
    factory.insertVertex({1,1,0});
    factory.insertVertex({0,0,1});
    factory.insertVertex({1,0,1});
    factory.insertVertex({0,1,1});
    factory.insertVertex({2,2,2});
    factory.insertElement(Dune::GeometryTypes::cube(3), {0,1,2,3,4,5,6,7});
    auto grid = factory.createGrid();
    grid->setRefinementType(Grid::RefinementType::COPY);
    refineLocalNearOrigin(*grid, 4);
    testSuite.subTest(checkOnGrid(*grid, "UGGrid<3> (nonaffine cubes)"));
  }

  return testSuite.exit();
}
catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
