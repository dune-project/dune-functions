// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>

using namespace Dune;

int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {4, 4};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisBuilder;

  auto basis = makeBasis(
    gridView,
    composite(
      pq<1>(),
      pq<1>(),
      pq<1>()
    ));

}
// Error handling
catch (Exception e)
{
  std::cout << e << std::endl;
}
