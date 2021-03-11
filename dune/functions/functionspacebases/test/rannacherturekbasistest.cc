// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/rannacherturekbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;
using namespace Dune::Functions;

template<int dim>
void testRannacherTurekBasis(TestSuite& test)
{
  ///////////////////////////
  /////  Simplex grids  /////
  ///////////////////////////

  // Check RannacherTurekBasis created 'manually'
  // Use grid with unknown-at-compile-time element type
  {
    using Grid = UGGrid<dim>;
    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    auto grid = (dim==2) ? GmshReader<Grid>::read(path + "curved2d.msh")
                         : GmshReader<Grid>::read(path + "telescope1storder.msh");

    Functions::RannacherTurekBasis<typename Grid::LeafGridView> basis(grid->leafGridView());
    test.subTest(checkBasis(basis, EnableCenterContinuityCheck()));
  }

  // Check RannacherTurekBasis created using basis builder mechanism
  // Use grid with unknown-at-compile-time element type
  {
    using Grid = UGGrid<dim>;
    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    auto grid = (dim==2) ? GmshReader<Grid>::read(path + "curved2d.msh")
                         : GmshReader<Grid>::read(path + "telescope1storder.msh");

    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid->leafGridView(), rannacherTurek());
    test.subTest(checkBasis(basis, EnableCenterContinuityCheck()));
  }


  ///////////////////////////
  /////    Cube grids   /////
  ///////////////////////////

  // Check RannacherTurekBasis created 'manually'
  // Use grid with known-at-compile-time element type
  {
    using Grid = YaspGrid<dim>;
    Dune::FieldVector<double,dim> one(1);
    std::array<int,dim> elems;
    elems.fill(5);

    Grid grid(one, elems);

    Functions::RannacherTurekBasis<typename Grid::LeafGridView> basis(grid.leafGridView());
    test.subTest(checkBasis(basis, EnableCenterContinuityCheck()));
  }

  // Check RannacherTurekBasis created using basis builder mechanism
  // Use grid with known-at-compile-time element type
  {
    using Grid = YaspGrid<dim>;
    Dune::FieldVector<double,dim> one(1);
    std::array<int,dim> elems;
    elems.fill(5);
    Grid grid(one, elems);

    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid.leafGridView(), rannacherTurek());
    test.subTest(checkBasis(basis, EnableCenterContinuityCheck()));
  }


  ///////////////////////////
  /////   Hybrid grids  /////
  ///////////////////////////

  // only for dim = 2
  // hybrid-testgrid-3d.msh contains pyramids and prisms, which are not implemented
  if (dim ==2)
  {
    // Check RannacherTurekBasis created 'manually'
    // Use grid with unknown-at-compile-time element type
    {
      using Grid = UGGrid<dim>;
      const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
      auto grid = (dim==2) ? GmshReader<Grid>::read(path + "hybrid-testgrid-2d.msh")
                         : GmshReader<Grid>::read(path + "hybrid-testgrid-3d.msh");

      Functions::RannacherTurekBasis<typename Grid::LeafGridView> basis(grid->leafGridView());
      test.subTest(checkBasis(basis, EnableCenterContinuityCheck()));
    }

    // Check RannacherTurekBasis created using basis builder mechanism
    // Use grid with unknown-at-compile-time element type
    {
      using Grid = UGGrid<dim>;
      const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
      auto grid = (dim==2) ? GmshReader<Grid>::read(path + "hybrid-testgrid-2d.msh")
                         : GmshReader<Grid>::read(path + "hybrid-testgrid-3d.msh");

      using namespace Functions::BasisFactory;
      auto basis = makeBasis(grid->leafGridView(), rannacherTurek());
      test.subTest(checkBasis(basis, EnableCenterContinuityCheck()));
    }
  }

}

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;


  testRannacherTurekBasis<2>(test);
  testRannacherTurekBasis<3>(test);


  return test.exit();

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
