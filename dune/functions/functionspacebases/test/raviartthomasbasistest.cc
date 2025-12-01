// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/functions/functionspacebases/raviartthomasbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template<int k, class GridFactory>
void testRaviartThomasBasis(TestSuite& test, const GridFactory& factory)
{
  auto grid = factory();
  auto gridView = grid->leafGridView();

  std::cout<<"  Testing order: "<< k <<std::endl;

  // Check RaviartThomasBasis created 'manually'
  Functions::RaviartThomasBasis<decltype(gridView),k> basis1(gridView);
  test.subTest(checkBasis(basis1, EnableNormalContinuityCheck()));

  // Check RaviartThomasBasis created using basis builder mechanism
  using namespace Functions::BasisFactory;
  auto basis2 = makeBasis(gridView, raviartThomas<k>());
  test.subTest(checkBasis(basis2, EnableNormalContinuityCheck()));

  // Now modify the grid, and check again.
  const auto firstEntity = gridView.template begin<0>();
  grid->mark(1, *firstEntity);
  grid->adapt();

  auto modifiedGridView = grid->leafGridView();

  // Check the RaviartThomasBasis that was created 'manually'
  basis1.update(modifiedGridView);
  test.subTest(checkBasis(basis1, EnableNormalContinuityCheck()));

  // Check the RaviartThomasBasis that was created using the basis builder mechanism
  basis2.update(modifiedGridView);
  test.subTest(checkBasis(basis2, EnableNormalContinuityCheck()));
}


int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

  TestSuite test;

  // Test with grid that only supports cube elements
  // (Grids with only a single element type receive special treatment by RaviartThomasBasis)
  std::cout<<"Testing RaviartThomasBasis in 2D with cube grids\n";

  auto quadGridFactory = []() {
    return std::make_unique<YaspGrid<2> >(FieldVector<double,2>{1.0, 1.0}, std::array<int,2>{5,5});
  };

  testRaviartThomasBasis<0>(test, quadGridFactory);
  testRaviartThomasBasis<1>(test, quadGridFactory);
  testRaviartThomasBasis<2>(test, quadGridFactory);

  std::cout<<"Testing RaviartThomasBasis in 3D with cube grids\n";

  auto hexaGridFactory = []() {
    return std::make_unique<YaspGrid<3> >(FieldVector<double,3>{1.0, 1.0, 1.0}, std::array<int,3>{4,4,4});
  };

  testRaviartThomasBasis<0>(test, hexaGridFactory);
  testRaviartThomasBasis<1>(test, hexaGridFactory);

  // Test with pure simplex grid
  // (Unfortunately there is no grid implementation available that only supports simplices.)
  std::cout<<"Testing RaviartThomasBasis in 2D with simplex grid\n";

  auto triangleGridFactory = [&path]() {
    return GmshReader<UGGrid<2> >::read(path + "curved2d.msh");
  };

  testRaviartThomasBasis<0>(test, triangleGridFactory);
  testRaviartThomasBasis<1>(test, triangleGridFactory);

  std::cout<<"Testing RaviartThomasBasis in 3D with simplex grid\n";

  auto tetraGridFactory = [&path]() {
    return GmshReader<UGGrid<3> >::read(path + "telescope.msh");
  };

  testRaviartThomasBasis<0>(test, tetraGridFactory);

  // Test with mixed-element 2d grid
  std::cout<<"Testing RaviartThomasBasis in 2D with mixed-element grid\n";

  auto mixed2dGridFactory = [&path]() {
    return GmshReader<UGGrid<2> >::read(path + "hybrid-testgrid-2d.msh");
  };

  testRaviartThomasBasis<0>(test, mixed2dGridFactory);
  //testRaviartThomasBasis<1>(test, mixed2dGridFactory);

  // Test with mixed-element 3d grid
  std::cout<<"Testing RaviartThomasBasis in 3D with mixed-element grid\n";

  auto mixed3dGridFactory = [&path]() {
    return GmshReader<UGGrid<3> >::read(path + "hybrid-testgrid-3d.msh");
  };

  testRaviartThomasBasis<0>(test, mixed3dGridFactory);

  return test.exit();
}
