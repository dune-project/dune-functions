// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/periodicbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
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
  using Grid = YaspGrid<dim>;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{4, 4}};
  Grid grid(l,elements);

  using GridView = Grid::LeafGridView;
  GridView gridView = grid.leafGridView();

  std::cout << "Host grid has " << gridView.size(dim) << " vertices." << std::endl;

  //////////////////////////////////////////////////
  //   Infrastructure for handling periodicity
  //////////////////////////////////////////////////

  // Check whether two points are equal on R/Z x R/Z
  auto equivalent = [](const FieldVector<double,2>& x, const FieldVector<double,2>& y)
  {
    return ( (FloatCmp::eq(x[0],y[0]) or FloatCmp::eq(x[0]+1,y[0]) or FloatCmp::eq(x[0]-1,y[0]))
           and (FloatCmp::eq(x[1],y[1]) or FloatCmp::eq(x[1]+1,y[1]) or FloatCmp::eq(x[1]-1,y[1])) );
  };

  /////////////////////////////////////////////////////////
  //   Test a PeriodicBasis all by itself
  /////////////////////////////////////////////////////////

  using namespace Functions::BasisFactory;
  using namespace Functions::BasisFactory::Experimental;

  {
    PeriodicIndexSet periodicIndices;

    // Don't do the following in real life: It has quadratic run-time in the number of vertices.
    for (const auto& v1 : vertices(gridView))
      for (const auto& v2 : vertices(gridView))
        if (equivalent(v1.geometry().corner(0), v2.geometry().corner(0)))
          periodicIndices.unifyIndexPair({gridView.indexSet().index(v1)}, {gridView.indexSet().index(v2)});

    auto basis = makeBasis(gridView, lagrange<1>());
    {
      auto periodicBasis = makeBasis(gridView, periodic(lagrange<1>(), periodicIndices));
      std::cout << "Solitary periodic basis has " << periodicBasis.dimension() << " degrees of freedom." << std::endl;
      test.subTest(checkBasis(periodicBasis, EnableContinuityCheck()));
    }
    {
      auto periodicBasis = makeBasis(gridView, periodic(basis.preBasis(), periodicIndices));
      std::cout << "Solitary periodic basis has " << periodicBasis.dimension() << " degrees of freedom." << std::endl;
      test.subTest(checkBasis(periodicBasis, EnableContinuityCheck()));
    }
    {
      auto periodicBasis = makeBasis(gridView, periodic(basis, periodicIndices));
      std::cout << "Solitary periodic basis has " << periodicBasis.dimension() << " degrees of freedom." << std::endl;
      test.subTest(checkBasis(periodicBasis, EnableContinuityCheck()));
    }
  }

  /////////////////////////////////////////////////////////
  //   Test a PeriodicBasis in a power basis
  /////////////////////////////////////////////////////////

  {
    PeriodicIndexSet periodicIndices;

    // Don't do the following in real life: It has quadratic run-time in the number of vertices.
    for (const auto& v1 : vertices(gridView))
      for (const auto& v2 : vertices(gridView))
        if (equivalent(v1.geometry().corner(0), v2.geometry().corner(0)))
          periodicIndices.unifyIndexPair({gridView.indexSet().index(v1)}, {gridView.indexSet().index(v2)});

    // Check whether power<periodic> does at least compile
    auto periodicBasis = makeBasis(
      gridView,
      power<dim>(
        periodic(lagrange<1>(), periodicIndices),
        blockedInterleaved()
    ));

    std::cout << "power<periodic> basis has " << periodicBasis.dimension() << " degrees of freedom" << std::endl;

    test.subTest(checkBasis(periodicBasis, EnableContinuityCheck()));
  }
}
catch (Exception& e)
{
  std::cout << e.what() << std::endl;
}
