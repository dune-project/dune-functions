// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stringutility.hh>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>

#include <dune/functions/functionspacebases/test/cubichermitebasis.hh>
#include <dune/functions/functionspacebases/test/interpolatetest.hh>
#include <dune/common/timer.hh>



using namespace Dune::Functions;



int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite suite;

  using namespace Dune::Functions::BasisFactory;

  {
    auto gridPtr = Dune::StructuredGridFactory<Dune::OneDGrid>::createCubeGrid({0}, {1}, {10});

    auto gridView = gridPtr->leafGridView();
    auto basis = makeBasis(gridView, cubicHermite());

//    suite.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
    suite.subTest(checkBasis(basis, EnableContinuityCheck()));

    {
      using Domain = Dune::FieldVector<double,1>;
      using Range = Dune::FieldVector<double,1>;
      using Jacobian = Dune::FieldMatrix<double,1,1>;

      using SignatureTag = Dune::Functions::SignatureTag<Range(Domain)>;
      auto f = [](const Domain& x) -> Range {
        return x*x;
      };
      auto df = [](const Domain& x) -> Jacobian {
        return Jacobian{2.0*x};
      };
      auto ff = Dune::Functions::makeDifferentiableFunctionFromCallables(SignatureTag(), f, df);
      auto coefficients = std::vector<double>();
      Dune::Functions::interpolate(basis, coefficients, ff);

      suite.subTest(checkInterpolateConsistency<Range>(basis, coefficients));
    }
  }
  {
    auto factory = Dune::GridFactory<Dune::UGGrid<2>>{};
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({2,0});
    factory.insertVertex({1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {0,1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {1,2,3});
    auto gridPtr = factory.createGrid();
    auto& grid = *gridPtr;
#if NDEBUG
    grid.globalRefine(7);
#else
    grid.globalRefine(3);
#endif

    using Domain = Dune::FieldVector<double,2>;
    using Range = Dune::FieldVector<double,1>;
    using Jacobian = Dune::FieldMatrix<double,1,2>;
    using SignatureTag = Dune::Functions::SignatureTag<Range(Domain)>;

    auto f = [](const auto& x) {
      return x[0]*x[0] + x[1]*x[1];
    };
    auto df = [](const auto& x) {
      return Jacobian({{2*x[0],2*x[1]}});
    };
    auto ff = Dune::Functions::makeDifferentiableFunctionFromCallables(SignatureTag(), f, df);


    {
      auto basis = makeBasis(grid.leafGridView(), cubicHermite());
//      suite.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
      suite.subTest(checkBasis(basis, EnableContinuityCheck()));

      auto coefficients = std::vector<double>();
      Dune::Functions::interpolate(basis, coefficients, ff);
      suite.subTest(checkInterpolateConsistency<Range>(basis, coefficients));
    }

    {
      auto basis = makeBasis(grid.leafGridView(), reducedCubicHermite());
//      suite.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
      suite.subTest(checkBasis(basis, EnableContinuityCheck()));

      auto coefficients = std::vector<double>();
      Dune::Functions::interpolate(basis, coefficients, ff);
      suite.subTest(checkInterpolateConsistency<Range>(basis, coefficients));
    }

  }

  return suite.exit();
}
