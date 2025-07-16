// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <array>
#include <cmath>
#include <numbers>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/functions/common/subdomain.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>
#include <dune/functions/functionspacebases/restrictedbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/rannacherturekbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/composedgridfunction.hh>
#include <dune/functions/backends/istlvectorfactory.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  // Create a Grid for testing
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

  grid->globalRefine(3);

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Functions::Experimental::BasisFactory;
  using namespace Dune::Indices;
  using Dune::Functions::Experimental::SubDomain;
  using Dune::Functions::interpolate;

  // Indicator functions for four square subdomains covering the four quadrants
  // and one overlapping square subdomain in the center
  auto subDomainAIndicator = [](auto x) {
    return (x[0] <= 0.5) and (x[1] <= 0.5);
  };

  auto subDomainBIndicator = [](auto x) {
    return (x[0] > 0.5) and (x[1] <= 0.5);
  };

  auto subDomainCIndicator = [](auto x) {
    return (x[0] <= 0.5) and (x[1] > 0.5);
  };

  auto subDomainDIndicator = [](auto x) {
    return (x[0] > 0.5) and (x[1] > 0.5);
  };

  auto subDomainEIndicator = [&](auto x) {
    return  (std::fabs(x[0] - 0.5) <= 0.25) and (std::fabs(x[1] - 0.5) <= 0.25);
  };

  // Create corresponding SubDomain objects
  auto subDomainA = SubDomain(gridView);
  for(auto&& element : elements(gridView))
    if (subDomainAIndicator(element.geometry().center()))
      subDomainA.insertElement(element);

  auto subDomainB = SubDomain(gridView);
  for(auto&& element : elements(gridView))
    if (subDomainBIndicator(element.geometry().center()))
      subDomainB.insertElement(element);

  auto subDomainC = SubDomain(gridView);
  for(auto&& element : elements(gridView))
    if (subDomainCIndicator(element.geometry().center()))
      subDomainC.insertElement(element);

  auto subDomainD = SubDomain(gridView);
  for(auto&& element : elements(gridView))
    if (subDomainDIndicator(element.geometry().center()))
      subDomainD.insertElement(element);

  auto subDomainE = SubDomain(gridView);
  for(auto&& element : elements(gridView))
    if (subDomainEIndicator(element.geometry().center()))
      subDomainE.insertElement(element);

  // Create a basis
  auto basis = makeBasis(gridView,
      composite(
        restrict(lagrange<3>(), subDomainA),
        restrict(
          composite(
            power<dim>(lagrange<2>(),flatInterleaved()),
            lagrange<1>(),
          flatLexicographic()),
          subDomainB),
        restrict(rannacherTurek(), subDomainC),
        restrict(lagrangeDG<5>(), subDomainD),
        lagrange<2>(),
        restrict(lagrange<5>(), subDomainE),
        blockedLexicographic()
      )
    );

  // Run basis check
  test.subTest(checkBasis(basis));

  // Now test grid functions on the composed basis:

  // We interpolate the sinsin function in the four quadrants
  // and once globally. Furthermore we interpolate a quartic bubble
  // on the center subdomain.
  //
  auto sinsin = [](auto x) {
    return std::sin(x[0]*2*std::numbers::pi)*std::sin(x[1]*2*std::numbers::pi);
  };

  auto bubble = [&] (auto x) {
    return 256*(x[0]-0.25)*(0.75-x[0])*(x[1]-0.25)*(0.75-x[1]);
  };

  auto c = Dune::Functions::makeISTLVector<double>(basis.preBasis().containerDescriptor());
  auto c_backend = Dune::Functions::istlVectorBackend(c);

  interpolate(subspaceBasis(basis, _0), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _1, _1), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _2), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _3), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _4), c_backend, sinsin);
  interpolate(subspaceBasis(basis, _5), c_backend, bubble);

  // Now we add the functions from all subdomains and subtract
  // the global sinsin interpolation. The result should be
  // (up to mismatch of the FE spaces) the bubble again.
  using Range = Dune::TupleVector<double, Dune::TupleVector<std::array<double,2>, double>, double, double, double, double, double>;
  auto c_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis, c_backend);
  auto bubble_gf = Dune::Functions::makeComposedGridFunction([](auto y) {
    return y[_0] + y[_1][_1] + y[_2] + y[_3] - y[_4] + y[_5];
  }, c_gf);

  // Finally we write all grid functions to vtk
  auto c0_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _0), c_backend);
  auto c1_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _1,_1), c_backend);
  auto c2_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _2), c_backend);
  auto c3_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _3), c_backend);
  auto c4_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _4), c_backend);
  auto c5_gf = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _5), c_backend);

  Dune::SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, Dune::refinementLevels(5));
  vtkWriter.addVertexData(c0_gf, Dune::VTK::FieldInfo("c0", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c1_gf, Dune::VTK::FieldInfo("c1", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c2_gf, Dune::VTK::FieldInfo("c2", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c3_gf, Dune::VTK::FieldInfo("c3", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c4_gf, Dune::VTK::FieldInfo("c4", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(c5_gf, Dune::VTK::FieldInfo("c5", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(bubble_gf, Dune::VTK::FieldInfo("bubble", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("restrictedprebasistest");

  return test.exit();
}
