// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/functions/functionspacebases/raviartthomasbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template<int k, class GridView>
void testRaviartThomasBasis(TestSuite& test, const GridView& gridView)
{
  std::cout<<"  Testing order: "<< k <<std::endl;

  // Check RaviartThomasBasis created 'manually'
  {
    Functions::RaviartThomasBasis<GridView,k> basis(gridView);
    test.subTest(checkBasis(basis, EnableNormalContinuityCheck()));
  }

  // Check RaviartThomasBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(gridView, raviartThomas<k>());
    test.subTest(checkBasis(basis, EnableNormalContinuityCheck()));
  }
}


int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  // Test with grid that only supports cube elements
  // (Grids with only a single element type receive special treatment by RaviartThomasBasis)
  std::cout<<"Testing RaviartThomasBasis in 2D with cube grids\n";
  YaspGrid<2> quadGrid({1.0, 1.0}, {5,5});
  auto quadGridView = quadGrid.leafGridView();

  testRaviartThomasBasis<0>(test, quadGridView);
  testRaviartThomasBasis<1>(test, quadGridView);
  testRaviartThomasBasis<2>(test, quadGridView);

  std::cout<<"Testing RaviartThomasBasis in 3D with cube grids\n";
  YaspGrid<3> hexaGrid({1.0, 1.0, 1.0}, {4,4,4});
  auto hexaGridView = hexaGrid.leafGridView();
  testRaviartThomasBasis<0>(test, hexaGridView);
  testRaviartThomasBasis<1>(test, hexaGridView);

  // Test with pure simplex grid
  // (Unfortunately there is no grid implementation available that only supports simplices.)
  std::cout<<"Testing RaviartThomasBasis in 2D with simplex grid\n";
  using Mixed2dGrid = UGGrid<2>;
  const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
  auto triangleGrid = GmshReader<Mixed2dGrid>::read(path + "curved2d.msh");
  auto triangleGridView = triangleGrid->leafGridView();
  testRaviartThomasBasis<0>(test, triangleGridView);
  testRaviartThomasBasis<1>(test, triangleGridView);

  std::cout<<"Testing RaviartThomasBasis in 3D with simplex grid\n";
  using Mixed3dGrid = UGGrid<3>;
  auto tetraGrid = GmshReader<Mixed3dGrid>::read(path + "telescope.msh");
  auto tetraGridView = tetraGrid->leafGridView();
  testRaviartThomasBasis<0>(test, tetraGridView);

  // Test with mixed-element 2d grid
  std::cout<<"Testing RaviartThomasBasis in 2D with mixed-element grid\n";
  auto mixed2dGrid = GmshReader<Mixed2dGrid>::read(path + "hybrid-testgrid-2d.msh");
  auto mixed2dGridView = mixed2dGrid->leafGridView();
  testRaviartThomasBasis<0>(test, mixed2dGridView);
  //testRaviartThomasBasis<1>(test, mixed2dGridView);

  return test.exit();
}
