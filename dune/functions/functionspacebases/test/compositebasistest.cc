// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/bitsetvector.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{4, 4}};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisFactory;

  auto basis = makeBasis(
    gridView,
    composite(
      lagrange<1>(),
      lagrange<1>(),
      lagrange<1>()
    ));

  test.subTest(checkBasis(basis, EnableContinuityCheck()));

  return test.exit();
}
// Error handling
catch (Exception& e)
{
  std::cout << e.what() << std::endl;
}
