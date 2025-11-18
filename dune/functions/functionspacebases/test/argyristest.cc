// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/timer.hh>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/argyrisbasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/test/enabledifferentiabilitycheck.hh>

using namespace Dune;
using namespace Dune::Functions;

int main(int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test("ArgyrisBasis");

  using Grid = UGGrid<2>;
  // using Grid = YaspGrid<2>;
  auto grid = StructuredGridFactory<Grid>::createSimplexGrid({0., 0.}, {1., 1.}, {3, 3});

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, argyris());

  test.subTest(checkBasis(basis, EnableContinuityCheck(),
                          EnableDifferentiabilityCheck(),
                          CheckLocalFiniteElementFlag<2>()));
  return test.exit();
}
