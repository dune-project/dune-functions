// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/nedelecbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template <int dim, int kind, int order>
void testNedelecBasis(TestSuite& test)
{
  ///////////////////////////
  /////  Simplex grids  /////
  ///////////////////////////

  // Check NedelecBasis created 'manually'
  // Use grid with unknown-at-compile-time element type
  {
    using Grid = UGGrid<dim>;
    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    auto grid = (dim==2) ? GmshReader<Grid>::read(path + "curved2d.msh")
                         : GmshReader<Grid>::read(path + "telescope1storder.msh");

    Functions::NedelecBasis<typename Grid::LeafGridView, kind, order, double> basis(grid->leafGridView());
    test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
  }

  // Check NedelecBasis created using basis builder mechanism
  // Use grid with unknown-at-compile-time element type
  {
    using Grid = UGGrid<dim>;
    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    auto grid = (dim==2) ? GmshReader<Grid>::read(path + "curved2d.msh")
                         : GmshReader<Grid>::read(path + "telescope1storder.msh");

    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid->leafGridView(), nedelec<kind,order, double>());
    test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
  }

  ///////////////////////////
  /////    Cube grids   /////
  ///////////////////////////

  // Check NedelecBasis created 'manually'
  // Use grid with known-at-compile-time element type
  {
    using Grid = YaspGrid<dim>;
    Dune::FieldVector<double,dim> one(1);
    std::array<int,dim> elems;
    elems.fill(5);

    Grid grid(one, elems);

    Functions::NedelecBasis<typename Grid::LeafGridView, kind, order, double> basis(grid.leafGridView());
    test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
  }

  // Check NedelecBasis created using basis builder mechanism
  // Use grid with known-at-compile-time element type
  {
    using Grid = YaspGrid<dim>;
    Dune::FieldVector<double,dim> one(1);
    std::array<int,dim> elems;
    elems.fill(5);
    Grid grid(one, elems);

    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid.leafGridView(), nedelec<kind,order, double>());
    test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
  }

  ///////////////////////////
  /////   Hybrid grids  /////
  ///////////////////////////

  // only for dim = 2
  // hybrid-testgrid-3d.msh contains pyramids and prisms, which are not implemented
  if (dim ==2)
  {
    // Check NedelecBasis created 'manually'
    // Use grid with unknown-at-compile-time element type
    {
      using Grid = UGGrid<dim>;
      const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
      auto grid = (dim==2) ? GmshReader<Grid>::read(path + "hybrid-testgrid-2d.msh")
                         : GmshReader<Grid>::read(path + "hybrid-testgrid-3d.msh");

      Functions::NedelecBasis<typename Grid::LeafGridView, kind, order, double> basis(grid->leafGridView());
      test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
    }

    // Check NedelecBasis created using basis builder mechanism
    // Use grid with unknown-at-compile-time element type
    {
      using Grid = UGGrid<dim>;
      const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
      auto grid = (dim==2) ? GmshReader<Grid>::read(path + "hybrid-testgrid-2d.msh")
                         : GmshReader<Grid>::read(path + "hybrid-testgrid-3d.msh");

      using namespace Functions::BasisFactory;
      auto basis = makeBasis(grid->leafGridView(), nedelec<kind,order, double>());
      test.subTest(checkBasis(basis, EnableTangentialContinuityCheck()));
    }
  }
}

int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  testNedelecBasis<2, 1, 1>(test);
  testNedelecBasis<3, 1, 1>(test);

  return test.exit();
}
