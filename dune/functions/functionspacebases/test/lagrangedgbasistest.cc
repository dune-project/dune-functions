// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;



  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);



  // check LagrangeDGBasis created 'manually'
  {
    typedef GridType::LeafGridView GridView;
    const GridView& gridView = grid.leafGridView();
    LagrangeDGBasis<GridView,2> basis(gridView);
    test.subTest(checkBasis(basis));
  }



  // check LagrangeDGBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid.leafGridView(), lagrangeDG<2>());
    test.subTest(checkBasis(basis));
  }



  return test.exit();
}
