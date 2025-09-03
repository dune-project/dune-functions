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
#include <dune/grid/onedgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/cubichermitebasis.hh>

using namespace Dune;
using namespace Dune::Functions;


int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test_1d("1d"), test_2d("2d"), test_3d("3d");

  using namespace Dune::Functions::BasisFactory;

  { // 1d
    std::cout << "CubicHermite test in 1d" << std::endl;
    auto grid = StructuredGridFactory<OneDGrid>::createSimplexGrid({0.}, {1.}, {10});

    auto gridView = grid->levelGridView(0);

    {
      std::cout << "Grid has " << gridView.size(0) << " elements and " << gridView.size(1)
                << " facets and " << gridView.size(2) << " vertices" << std::endl;

      auto basis = makeBasis(gridView, cubicHermite());
      std::cout << "Basis has " << basis.size() << " dofs" << std::endl;

      test_1d.subTest(checkBasis(basis, EnableContinuityCheck(), CheckLocalFiniteElementFlag<2>{}, EnableDifferentiabilityCheck(),
                                 EnableVertexDifferentiabilityCheck()));
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
    std::cout << "Grid has " << gridView.size(0) << " elements and " << gridView.size(1)
              << " facets and " << gridView.size(2) << " vertices" << std::endl;

    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, cubicHermite());
      std::cout << "Basis has " << basis.size() << " dofs" << std::endl;

      test_2d.subTest(checkBasis(basis, CheckLocalFiniteElementFlag<1>{}, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck()));
    }
  }

  { // 2d  reduced
    std::cout << "reduced CubicHermite test in 2d" << std::endl;

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
    std::cout << "Grid has " << gridView.size(0) << " elements and " << gridView.size(1)
              << " facets and " << gridView.size(2) << " vertices" << std::endl;
    {
      using namespace Dune::Functions::BasisFactory;

      auto basis = makeBasis(gridView, reducedCubicHermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_2d.subTest(checkBasis(basis,CheckLocalFiniteElementFlag<1>{}, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck()));
    }
  }

  { // 3d
    std::cout << "CubicHermite test in 3d" << std::endl;

    auto grid = StructuredGridFactory<UGGrid<3>>::createSimplexGrid({0., 0., 0.}, {1., 1., 1.},
                                                                    {{3, 3, 3}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elements and " << gridView.size(1)
              << " facets and " << gridView.size(2) << " edges and " << gridView.size(3)
              << " vertices " << std::endl;
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, cubicHermite());
      std::cout << "Basis has " << basis.size() << " dofs" << std::endl;

      test_3d.subTest(checkBasis(basis, EnableContinuityCheck(),CheckLocalFiniteElementFlag<1>{},
                                 EnableVertexDifferentiabilityCheck()));
    }
  }

  return test_1d.exit() + test_2d.exit() + test_3d.exit();
}
