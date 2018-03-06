// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/raviartthomasbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;



  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);



  // check RaviartThomasBasis created 'manually'
  {
    typedef GridType::LeafGridView GridView;
    const GridView& gridView = grid.leafGridView();
    Functions::RaviartThomasBasis<GridView,0,GeometryType::BasicType::cube> basis(gridView);
    test.subTest(checkBasis(basis));
  }



  // check RaviartThomasBasis created using basis builder mechanism
  {
    using namespace Functions::BasisBuilder;
    auto basis = makeBasis(grid.leafGridView(), rt<0, GeometryType::BasicType::cube>());
    test.subTest(checkBasis(basis));
  }

  return test.exit();
}
