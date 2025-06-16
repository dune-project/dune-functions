// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>
#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/dynamicpowerbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>



// This wrapper will hide the derivative() of an existing GridFunction
// making it non-differentiable for testing.
template<class GlobalBase>
class NonDifferentiableGridFunction : protected GlobalBase
{
  using LocalBase = std::decay_t<decltype(localFunction(std::declval<GlobalBase>()))>;

  class LocalFunction : private LocalBase
  {
  public:
    LocalFunction(const LocalBase& base) : LocalBase(base) {}
    LocalFunction(LocalBase&& base) : LocalBase(base) {}
    using LocalBase::operator();
    using LocalBase::bind;
    using LocalBase::unbind;
    using LocalBase::bound;
    using LocalBase::localContext;
  };

public:
  NonDifferentiableGridFunction(GlobalBase&& base) : GlobalBase(base) {}
  using GlobalBase::operator();
  using GlobalBase::entitySet;
  friend LocalFunction localFunction(const NonDifferentiableGridFunction& f)
  {
    return LocalFunction(localFunction(static_cast<const GlobalBase&>(f)));
  }
};



template<class GridView>
Dune::TestSuite checkMakeBasis(const GridView& gridView)
{
  Dune::TestSuite test;

  using namespace Dune::Functions::BasisFactory;

  const int N = 10;
  const int M = 10;

  auto basis = makeBasis(gridView,
      power<N>(       // static power node
        power(        // dynamic power node
          composite(
            lagrange<3>(),  // lagrange basis with static order
            lagrange(1),    // lagrange basis with dynamic order
            flatLexicographic()),
          M,
          flatLexicographic()),
        blockedInterleaved())
      );

  test.subTest(checkBasis(basis, EnableContinuityCheck()));

  {
    [[maybe_unused]] auto& firstLagrangeFactor = basis.preBasis().subPreBasis().subPreBasis().subPreBasis(Dune::Indices::_0);
    [[maybe_unused]] auto& secondLagrangeFactor = basis.preBasis().subPreBasis().subPreBasis().template subPreBasis<1>();
  }

  using Vector = std::vector<Dune::FieldVector<double,N>>;

  Vector x;

  auto f = [](const auto& x){
    std::array<std::array<Dune::FieldVector<double,2>, M>, N> y;
    for(auto& yi : y)
      for(auto& yij : yi)
        yij = 1.0;
    return y;
  };

  Dune::Functions::interpolate(basis, x, f);

  for(const auto& xi : x)
    for(const auto& xij : xi)
      test.require(std::abs(xij - 1.0) < 1e-10)
        << "Coefficient of interpolated 1-function does not match";

  // Now check the same but provide f as non-differentiable GridFunction
  auto nonDiffF = NonDifferentiableGridFunction(Dune::Functions::makeAnalyticGridViewFunction(f, basis.gridView()));
  Dune::Functions::interpolate(basis, x, nonDiffF);

  for(const auto& xi : x)
    for(const auto& xij : xi)
      test.require(std::abs(xij - 1.0) < 1e-10)
        << "Coefficient of interpolated 1-function does not match";

  {
    auto basis = makeBasis(gridView,
        power<2>(
          composite(
            power<1>(power<1>(lagrange<1>(), blockedInterleaved()), blockedInterleaved()),
            power<2>(lagrange<1>(), blockedInterleaved()),
            power<3>(lagrange<1>(), blockedInterleaved()),
            blockedLexicographic()
          ),
          flatInterleaved()
        )
      );
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    auto basis = makeBasis(gridView,
        power<2>(
          composite(
            power<1>(power<1>(lagrange<1>(), blockedInterleaved()), blockedInterleaved()),
            power<2>(lagrange<1>(), blockedInterleaved()),
            power<3>(lagrange<1>(), blockedInterleaved()),
            blockedLexicographic()
          ),
          flatLexicographic()
        )
      );
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  return test;
}


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  // Check with YaspGrid
  {
    const int dim = 2;
    using Grid = Dune::YaspGrid<dim>;
    Dune::FieldVector<double,dim> l(1);
    std::array<int,dim> elements = {{10, 10}};
    Grid grid(l,elements);

    auto gridView = grid.leafGridView();

    test.subTest(checkMakeBasis(gridView));
  }

  // Check with mixed UGGrid
  {
    using Grid = Dune::UGGrid<2>;

    auto factory = Dune::GridFactory<Grid>();
    for(unsigned int k : Dune::range(9))
      factory.insertVertex({0.5*(k%3), 0.5*(k/3)});
    factory.insertElement(Dune::GeometryTypes::cube(2), {0, 1, 3, 4});
    factory.insertElement(Dune::GeometryTypes::cube(2), {1, 2, 4, 5});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {3, 4, 6});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {4, 7, 6});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {4, 5, 7});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {5, 8, 7});

    auto gridPtr = std::unique_ptr(factory.createGrid());
    auto& grid = *gridPtr;

    auto gridView = grid.leafGridView();

    test.subTest(checkMakeBasis(gridView));
  }

  return test.exit();
}
