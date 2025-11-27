// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/functions/functionspacebases/brezzidouglasmarinibasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;


template <int order, class GridFactory>
void testBrezziDouglasMariniBasis(TestSuite& test, const GridFactory& factory)
{
  auto grid = factory();
  auto gridView = grid->leafGridView();

  std::cout<<"  Testing order: "<< order <<std::endl;

  // Check BrezziDouglasMariniBasis created 'manually'
  Functions::BrezziDouglasMariniBasis<decltype(gridView), order> basis1(gridView);
  test.subTest(checkBasis(basis1, EnableNormalContinuityCheck()));

  // Check BrezziDouglasMariniBasis created using basis builder mechanism
  using namespace Functions::BasisFactory;
  auto basis2 = makeBasis(gridView, brezziDouglasMarini<order>());
  test.subTest(checkBasis(basis2, EnableNormalContinuityCheck()));
}


int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

  TestSuite test;

  std::cout<<"Testing NedelecBasis in 2D with simplex grid\n";
  auto triangleGridFactory = [&path]() {
    return GmshReader<UGGrid<2> >::read(path + "curved2d.msh");
  };

  testBrezziDouglasMariniBasis<1>(test, triangleGridFactory);
  // TODO: Enable this test!
  //testBrezziDouglasMariniBasis<2>(test, triangleGridFactory);

  // Test with grid that only supports cube elements
  std::cout<<"Testing BrezziDouglasMariniBasis in 2D with cube grid\n";
  auto quadGridFactory = []() {
    return std::make_unique<YaspGrid<2> >(FieldVector<double,2>{1.0, 1.0}, std::array<int,2>{5,5});
  };

  testBrezziDouglasMariniBasis<1>(test, quadGridFactory);
  // TODO: Enable this test!
  //testBrezziDouglasMariniBasis<2>(test, quadGridFactory);

#if 0  // TODO: Enable this test!
  std::cout<<"Testing BrezziDouglasMariniBasis in 2D with mixed-element grid\n";
  auto mixed2dGridFactory = [&path]() {
    return GmshReader<UGGrid<2> >::read(path + "hybrid-testgrid-2d.msh");
  };

  testBrezziDouglasMariniBasis<1>(test, mixed2dGridFactory);
#endif

#if 0  // TODO: Enable this test!
  std::cout<<"Testing BrezziDouglasMariniBasis in 3D with simplex grid\n";

  auto tetraGridFactory = [&path]() {
    return GmshReader<UGGrid<3> >::read(path + "telescope1storder.msh");
  };

  testBrezziDouglasMariniBasis<1>(test, tetraGridFactory);
#endif

#if 0  // TODO: Enable this test!
  Test with grid that only supports cube elements
  std::cout<<"Testing BrezziDouglasMariniBasis in 3D with cube grid\n";
  auto cubeGridFactory = []() {
    return std::make_unique<YaspGrid<3> >(FieldVector<double,3>{1.0, 1.0, 1.0}, std::array<int,3>{5,5,5});
  };

  testBrezziDouglasMariniBasis<1>(test, cubeGridFactory);
#endif

#if 0  // TODO: Enable this test!
  hybrid-testgrid-3d.msh contains pyramids and prisms, which are not implemented
  std::cout<<"Testing BrezziDouglasMariniBasis in 3D with mixed-element grid\n";
  auto mixed3dGridFactory = [&path]() {
    return GmshReader<UGGrid<3> >::read(path + "hybrid-testgrid-3d.msh");
  };

  testBrezziDouglasMariniBasis<1>(test, mixed3dGridFactory);
#endif

  return test.exit();
}
