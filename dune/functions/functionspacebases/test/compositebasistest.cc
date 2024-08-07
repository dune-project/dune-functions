// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/bitsetvector.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{4, 4}};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisFactory;

  {
    auto basis = makeBasis(
      gridView,
      composite(
        lagrange<1>(),
        lagrange<1>(),
        lagrange<1>()
      ));

    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    auto basis = Dune::Functions::DefaultGlobalBasis(
      gridView,
      composite(
        lagrange<1>(),
        lagrange<1>(),
        lagrange<1>()
      ));

    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    using namespace Functions;
    using PreBasis =
      PowerPreBasis<BlockedInterleaved,
        CompositePreBasis<BlockedLexicographic,
          PowerPreBasis<BlockedInterleaved,
            LagrangePreBasis<GridView,1>,
            2>,
          PowerPreBasis<BlockedInterleaved,
            LagrangePreBasis<GridView,2>,
            2>
        >,
      2>;
    auto basis = Dune::Functions::DefaultGlobalBasis<PreBasis>(grid.leafGridView());
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    using namespace Functions;
    using PreBasis = CompositePreBasis<BlockedLexicographic, LagrangePreBasis<GridView,1>>;
    auto basis = Dune::Functions::DefaultGlobalBasis<PreBasis>(grid.leafGridView());
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    using namespace Functions;
    using PreBasis = CompositePreBasis<BlockedLexicographic, LagrangePreBasis<GridView,1>>;
    const auto gridView = grid.leafGridView();
    auto basis = Dune::Functions::DefaultGlobalBasis<PreBasis>(gridView);
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    using namespace Functions;
    using PreBasis = CompositePreBasis<BlockedLexicographic, LagrangePreBasis<GridView,1>>;
    auto basis = Dune::Functions::DefaultGlobalBasis<PreBasis>(gridView);
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  return test.exit();
}
// Error handling
catch (Exception& e)
{
  std::cout << e.what() << std::endl;
}
