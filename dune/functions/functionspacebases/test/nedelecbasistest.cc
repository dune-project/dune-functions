// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/nedelecbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template <int dim, int kind, int order>
void testNedelecBasis(TestSuite& test)
{
  // Check NedelecBasis created 'manually'
  // Use grid with unknown-at-compile-time element type
  {
    using Grid = UGGrid<dim>;
    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    auto grid = GmshReader<Grid>::read(path + "curved2d.msh");

    Functions::NedelecBasis<typename Grid::LeafGridView, kind, order, double> basis(grid->leafGridView());
    test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
  }

  // Check NedelecBasis created using basis builder mechanism
  // Use grid with unknown-at-compile-time element type
  {
    using Grid = UGGrid<dim>;
    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    auto grid = GmshReader<Grid>::read(path + "curved2d.msh");

    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid->leafGridView(), nedelec<kind,order, double>());
    test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
  }
}

int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  testNedelecBasis<2, 1, 1>(test);

  return test.exit();
}
