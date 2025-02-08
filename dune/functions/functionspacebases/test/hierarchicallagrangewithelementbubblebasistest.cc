// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/hierarchicallagrangewithelementbubblebasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template <class GridView>
TestSuite checkDimension (const GridView& gridView)
{
  TestSuite test("Check dimension " + std::to_string(GridView::dimension));

  { // check HierarchicalBasis created 'manually'
    Functions::HierarchicalLagrangeWithElementBubbleBasis<GridView,1> basis1(gridView);
    test.subTest(checkBasis(basis1, EnableContinuityCheck()));
    Functions::HierarchicalLagrangeWithElementBubbleBasis<GridView,2> basis2(gridView);
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
  }

  { // check HierarchicalBasis created using basis builder mechanism
    using namespace Functions::BasisFactory;
    auto basis1 = makeBasis(gridView, hierarchicalLagrangeWithElementBubble<1>());
    test.subTest(checkBasis(basis1, EnableContinuityCheck()));
    auto basis2 = makeBasis(gridView, hierarchicalLagrangeWithElementBubble<2>());
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
  }

  { // check HierarchicalBasis created using basis builder mechanism
    using namespace Functions::BasisFactory;
    auto basis1 = makeBasis(gridView, hierarchicalP1B());
    test.subTest(checkBasis(basis1, EnableContinuityCheck()));
    auto basis2 = makeBasis(gridView, hierarchicalP2B());
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
  }

  return test;
}

int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  { // 1d
    using Factory = StructuredGridFactory<OneDGrid>;
    auto grid = Factory::createSimplexGrid({0.0}, {1.0}, {4u});
    test.subTest(checkDimension(grid->leafGridView()));
  }

  { // 2d
    using Factory = StructuredGridFactory<UGGrid<2>>;
    auto grid = Factory::createSimplexGrid({0.0,0.0}, {1.0,1.0}, {4u,4u});
    test.subTest(checkDimension(grid->leafGridView()));
  }

  { // 3d
    using Factory = StructuredGridFactory<UGGrid<3>>;
    auto grid = Factory::createSimplexGrid({0.0,0.0,0.0}, {1.0,1.0,1.0}, {4u,4u,4u});
    test.subTest(checkDimension(grid->leafGridView()));
  }

  return test.exit();
}
