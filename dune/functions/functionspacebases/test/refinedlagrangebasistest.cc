// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/albertagrid.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/refinedlagrangebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace Dune;
using namespace Dune::Functions;



int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;



  using namespace Dune::Functions::BasisFactory;

  {
    const int dim = 2;
    using Grid = Dune::AlbertaGrid<dim,dim>;

    std::unique_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid({0.0,0.0},{1.0,1.0},{1u,1u});

    grid->globalRefine(2);

    auto gridView = grid->leafGridView();
    auto basis0 = makeBasis(gridView, refinedLagrange<0>());
    auto basis1 = makeBasis(gridView, refinedLagrange<1>());

    test.subTest(checkBasis(basis0));
    test.subTest(checkBasis(basis1, EnableContinuityCheck()));

    std::vector<double> v0,v1;
    v0.resize(basis0.size(), 0);
    v1.resize(basis1.size(), 0);
    v0[2] = 1;
    v1[5] = 1;
    auto v0_f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis0,v0);
    auto v1_f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis1,v1);

    SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, Dune::refinementLevels(5));
    vtkWriter.addVertexData(v0_f, VTK::FieldInfo("refinedP0", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(v1_f, VTK::FieldInfo("refinedP1", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("debug");

  }

  {
    std::unique_ptr<OneDGrid> grid
      = StructuredGridFactory<OneDGrid>::createSimplexGrid({0}, {1}, {10});

    auto gridView = grid->levelGridView(0);

    {
      auto basis0 = makeBasis(gridView, refinedLagrange<0>());
      test.subTest(checkBasis(basis0));
      auto basis1 = makeBasis(gridView, refinedLagrange<1>());
      test.subTest(checkBasis(basis1, EnableContinuityCheck()));
    }

    {
      auto basis0 = makeBasis(gridView, refinedLagrange<0,float>());
      test.subTest(checkBasis(basis0));
      auto basis1 = makeBasis(gridView, refinedLagrange<1,float>());
      test.subTest(checkBasis(basis1, EnableContinuityCheck()));
    }

  }


  return test.exit();
}
