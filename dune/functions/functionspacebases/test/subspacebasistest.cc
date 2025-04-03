// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>




int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  auto grid2d = [](){
    const int dim = 2;
    using Grid = Dune::UGGrid<dim>;

    Dune::GridFactory<Grid> factory;
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({1,0});
    factory.insertVertex({1,1});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {0,1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {1,2,3});

    auto grid = factory.createGrid();
    grid->globalRefine(2);
    return grid;
  }();

  auto sameSubspace = [](const auto& sb1, const auto& sb2) {
    return (((void*)&sb1.rootBasis())==((void*)&sb2.rootBasis())) and (sb1.prefixPath()._data==sb2.prefixPath()._data);
  };

  using namespace Dune;
  using namespace Dune::Functions;
  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Indices;

  {
    auto basis = makeBasis(grid2d->leafGridView(), power<3>(lagrange<3>(), blockedInterleaved()));
    auto tp = Dune::TypeTree::treePath(_0);
    auto sb0_a = subspaceBasis(basis, tp);
    auto sb0_b = subspaceBasis(basis, _0);
    test.check(sameSubspace(sb0_a, sb0_b)) << "SubspaceBasis created from treepath and treepath indices are different";
  }

  {
    auto basis = makeBasis(grid2d->leafGridView(), composite(power<3>(lagrange<3>(), blockedInterleaved()), lagrange<1>(), blockedLexicographic()));
    auto sb0 = subspaceBasis(basis, _0);
    auto sb01_a = subspaceBasis(sb0, 1);
    auto sb01_b = subspaceBasis(basis, _0, 1);
    test.check(sameSubspace(sb01_a, sb01_b)) << "subspaceBasis(subspaceBasis(b,tp1),tp2) not equals subspaceBasis(b,(tp1,tp2))";
  }

  {
    auto basis = makeBasis(grid2d->leafGridView(), composite(power<3>(lagrange<3>(), blockedInterleaved()), lagrange<1>(), blockedLexicographic()));
    auto sb0 = SubspaceBasis(basis, Dune::TypeTree::treePath(_0));
    auto sb01_a = SubspaceBasis(sb0, Dune::TypeTree::treePath(1));
    auto sb01_b = SubspaceBasis(basis, Dune::TypeTree::treePath(_0, 1));
    test.check(sameSubspace(sb01_a, sb01_b)) << "SubspaceBasis(SubspaceBasis(b,tp1),tp2) not equals SubspaceBasis(b,(tp1,tp2))";
  }

  {
    auto basis = makeBasis(grid2d->leafGridView(), composite(power<3>(lagrange<3>(), blockedInterleaved()), lagrange<1>(), blockedLexicographic()));
    auto sb0 = SubspaceBasis(basis, Dune::TypeTree::treePath(_0));
    auto sb01_a = SubspaceBasis(sb0, Dune::TypeTree::treePath(1));
    auto sb01_b = subspaceBasis(basis, _0, 1);
    test.check(sameSubspace(sb01_a, sb01_b)) << "SubspaceBasis(SubspaceBasis(b,tp1),tp2) not equals subspaceBasis(b,(tp1,tp2))";
  }



  return test.exit();
}
catch (Dune::Exception& e)
{
  std::cout << e.what() << std::endl;
}
