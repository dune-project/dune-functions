// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/hierarchicallagrangebasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template <class Grid>
void testDim(TestSuite& test)
{
  static const int dim = Grid::dimension;
  using Factory = StructuredGridFactory<Grid>;
  FieldVector<double,dim> lower(0.0), upper(1.0);
  std::array<unsigned int,dim> elems = Dune::filledArray<dim,unsigned int>(1<<(5-dim));
  auto gridPtr = Factory::createSimplexGrid(lower, upper, elems);

  using GridView = typename Grid::LeafGridView;
  auto gridView = gridPtr->leafGridView();

  // check HierarchicalBasis created 'manually'
  {
    Functions::HierarchicalLagrangeBasis<GridView,1> basis1(gridView);
    test.subTest(checkBasis(basis1, EnableContinuityCheck()));

    Functions::HierarchicalLagrangeBasis<GridView,2> basis2(gridView);
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
  }

  // check HierarchicalBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis1 = makeBasis(gridView, hierarchicalLagrange<1>());
    test.subTest(checkBasis(basis1, EnableContinuityCheck()));

    auto basis2 = makeBasis(gridView, hierarchicalLagrange<2>());
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
  }
}

int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;
  testDim<OneDGrid>(test);
  testDim<UGGrid<2>>(test);
  testDim<UGGrid<3>>(test);

  return test.exit();
}
