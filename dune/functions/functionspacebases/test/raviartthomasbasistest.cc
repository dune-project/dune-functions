// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/raviartthomasbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template<int dim,int k>
void testRaviartThomasCube(TestSuite& test)
{
  FieldVector<double,dim> l(1);
  FieldVector<double,dim> ll(0);

  std::cout<<"  order: "<<k<<std::endl;

  typedef YaspGrid<dim> GridType;
  std::array<int,dim> elements{};
  elements.fill(10);

  GridType grid(l,elements);
  // check RaviartThomasBasis created 'manually'
  {
    Functions::RaviartThomasBasis<typename GridType::LeafGridView,k> basis(grid.leafGridView());
    test.subTest(checkBasis(basis, EnableNormalContinuityCheck()));
  }

  // check RaviartThomasBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid.leafGridView(), raviartThomas<k>());
    test.subTest(checkBasis(basis, EnableNormalContinuityCheck()));
  }

  // check RaviartThomasBasis on a grid without a compile-time-fixed element type
  {
    using Grid = UGGrid<dim>;
    std::array<unsigned int,dim> elements{};
    elements.fill(10);

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createCubeGrid(ll, l,elements);
    Functions::RaviartThomasBasis<typename Grid::LeafGridView,k> basis(grid->leafGridView());
    test.subTest(checkBasis(basis, EnableNormalContinuityCheck()));
  }
}

template<int dim,int k>
void testRaviartThomasSimplex(TestSuite& test)
{
  FieldVector<double,dim> l(1);
  FieldVector<double,dim> ll(0);

  std::cout<<"  order: "<<k<<std::endl;

  typedef UGGrid<dim> GridType;
  std::array<unsigned int,dim> elements{};
  elements.fill(10);
  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(ll, l, elements);

  // check RaviartThomasBasis created 'manually'
  {
    Functions::RaviartThomasBasis<typename GridType::LeafGridView,k> basis(grid->leafGridView());
    test.subTest(checkBasis(basis, EnableNormalContinuityCheck()));
  }

  // check RaviartThomasBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(grid->leafGridView(), raviartThomas<k>());
    test.subTest(checkBasis(basis, EnableNormalContinuityCheck()));
  }
}

int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  std::cout<<"Testing RaviartThomasBasis in 2D with cube grids\n";
  testRaviartThomasCube<2,0>(test);
  testRaviartThomasCube<2,1>(test);
  //cant be tested because there is no simplex RaviartThomas implementation for dim=2,k=2
  //testRaviartThomasCube<2,2>(test);

  std::cout<<"Testing RaviartThomasBasis in 3D with cube grids\n";
  testRaviartThomasCube<3,0>(test);
  //cant be tested because there is no simplex RaviartThomas implementation for dim=3,k=1
  //testRaviartThomasCube<3,1>(test);

  std::cout<<"Testing RaviartThomasBasis in 2D with simplex grids\n";
  testRaviartThomasSimplex<2,0>(test);
  // seems to be broken
  //testRaviartThomasSimplex<2,1>(test);

  std::cout<<"Testing RaviartThomasBasis in 3D with simplex grids\n";
  testRaviartThomasSimplex<3,0>(test);

  return test.exit();
}
