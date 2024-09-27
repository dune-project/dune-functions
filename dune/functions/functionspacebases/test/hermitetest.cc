// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/hermitebasis.hh>
#include <dune/functions/functionspacebases/test/cubichermitebasis.hh>

using namespace Dune;
using namespace Dune::Functions;

// Hack: Disable test that has not been merged to master so far,
// by replacing it with a dummy.
template<int i=0>
class CheckLocalFiniteElementFlag {};

template<class LFE>
auto benchmarkEvaluation(LFE const &lfe, int repeat = 1000)
{
  using RT = typename LFE::Traits::LocalBasisType::Traits::RangeType;
  using DT = typename LFE::Traits::LocalBasisType::Traits::DomainType;
  std::vector<RT> out;
  out.resize(lfe.size());
  Dune::Timer timer;
  for (auto i : Dune::range(repeat))
    lfe.localBasis().evaluateFunction(DT{}, out);
  return timer.elapsed();
}

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test_1d("1d"), test_2d("2d"), test_3d("3d");

  using namespace Dune::Functions::BasisFactory;
  bool benchmark = true;
  int repeat = 1000000;

  { // 1d
    std::cout << "Hermite test in 1d" << std::endl;
    auto grid = StructuredGridFactory<OneDGrid>::createSimplexGrid({0.}, {1.}, {10});

    auto gridView = grid->levelGridView(0);

    {
      std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
                << " facettes and " << gridView.size(2) << " vertices" << std::endl;

      auto basis = makeBasis(gridView, hermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_1d.subTest(checkBasis(basis, EnableContinuityCheck(), EnableDifferentiabilityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag<2>()));
      if (benchmark) {
        auto lv = basis.localView();
        for (auto e : elements(gridView)) {
          lv.bind(e);
          std::cout << repeat << " Evaluations took "
                    << benchmarkEvaluation(lv.tree().finiteElement(), repeat) << std::endl;
          break;
        }
      }
    }
  }

    { // 1d
    std::cout << "Hermite test in 1d, Carstens implementation" << std::endl;
    auto grid = StructuredGridFactory<OneDGrid>::createSimplexGrid({0.}, {1.}, {10});

    auto gridView = grid->levelGridView(0);

    {
      std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
                << " facettes and " << gridView.size(2) << " vertices" << std::endl;

      auto basis = makeBasis(gridView, cubicHermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_1d.subTest(checkBasis(basis, EnableContinuityCheck(), EnableDifferentiabilityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag<0>()));
      if (benchmark) {
        auto lv = basis.localView();
        for (auto e : elements(gridView)) {
          lv.bind(e);
          std::cout << repeat << " Evaluations took "
                    << benchmarkEvaluation(lv.tree().finiteElement(), repeat) << std::endl;
          break;
        }
      }
    }
  }

  { // 2d
    std::cout << "Hermite test in 2d" << std::endl;

    auto factory = Dune::GridFactory<Dune::UGGrid<2>>{};
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({2,0});
    factory.insertVertex({1,4});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {0,1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {1,2,3});
    auto grid = factory.createGrid();
    grid->globalRefine(2);

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    // using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, hermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_2d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag<2>()));
      if (benchmark) {
        auto lv = basis.localView();
        for (auto e : elements(gridView)) {
          lv.bind(e);
          std::cout << repeat << " Evaluations took "
                    << benchmarkEvaluation(lv.tree().finiteElement(), repeat) << std::endl;
          break;
        }
      }
    }
  }

  { // 2d  reduced
    std::cout << "reduced Hermite test in 2d" << std::endl;

    auto factory = Dune::GridFactory<Dune::UGGrid<2>>{};
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({2,0});
    factory.insertVertex({1,4});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {0,1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {1,2,3});
    auto grid = factory.createGrid();
    grid->globalRefine(2);

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    {
      using namespace Dune::Functions::BasisFactory;

      auto basis = makeBasis(gridView, reducedHermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_2d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag<2>()));
      if (benchmark) {
        auto lv = basis.localView();
        for (auto e : elements(gridView)) {
          lv.bind(e);
          std::cout << repeat << " Evaluations took "
                    << benchmarkEvaluation(lv.tree().finiteElement(), repeat) << std::endl;
          break;
        }
      }
    }
  }

  { // 2d  reduced
    std::cout << "reduced Hermite test in 2d, Carstens Implementation" << std::endl;

    auto grid = StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0., 0.}, {1., 1.}, {{2, 2}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    {
      using namespace Dune::Functions::BasisFactory;

      auto basis = makeBasis(gridView, reducedCubicHermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_2d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag<0>()));
      if (benchmark) {
        auto lv = basis.localView();
        for (auto e : elements(gridView)) {
          lv.bind(e);
          std::cout << repeat << " Evaluations took "
                    << benchmarkEvaluation(lv.tree().finiteElement(), repeat)
                    << std::endl;
          break;
        }
      }
    }
  }

  { // 3d
    std::cout << "Hermite test in 3d" << std::endl;

    auto grid = StructuredGridFactory<UGGrid<3>>::createSimplexGrid({0., 0., 0.}, {1., 1., 1.},
                                                                    {{3, 3, 3}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " edges and " << gridView.size(3)
              << " vertices " << std::endl;
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, hermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_3d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag<2>()));
      if (benchmark) {
        auto lv = basis.localView();
        for (auto e : elements(gridView)) {
          lv.bind(e);
          std::cout << repeat << " Evaluations took "
                    << benchmarkEvaluation(lv.tree().finiteElement(), repeat) << std::endl;
          break;
        }
      }
    }
  }

  return test_1d.exit() + test_2d.exit() + test_3d.exit();
}
