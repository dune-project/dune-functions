// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/dynamicdgbasis.hh>

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
  typedef GridType::LeafGridView GridView;

  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);
  const GridView& gridView = grid.leafGridView();

  using Basis = Dune::Functions::DynamicQkDGBasis<GridView>;
  auto mapper = Basis::PreBasis::mapper(gridView);
  auto degrees = std::vector<int>(mapper.size());
  // set some arbitrary orders to demonstrate
  // the concept:
  for(size_t i = 0; i < degrees.size(); i++) {
    degrees[i] = i%2 +1;
  }

  // check DynamicDGBasis created 'manually'
  {
    auto basis = Basis{gridView, degrees};

    test.subTest(checkBasis(basis));
  }



  // check DynamicDGBasis created using basis builder mechanism
  {
    using namespace Functions::BasisBuilder;
    auto basis = makeBasis(grid.leafGridView(), dynamicDG(degrees));
    test.subTest(checkBasis(basis));
  }

  return test.exit();
}
