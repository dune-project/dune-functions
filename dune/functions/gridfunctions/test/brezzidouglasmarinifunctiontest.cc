// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/functions/functionspacebases/brezzidouglasmarinibasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>
#include <dune/functions/gridfunctions/test/continuitychecks.hh>


using namespace Dune;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
  bool passed = true;

  // Generate grid for testing
  const int dim = 2;
  using Grid = UGGrid<dim>;

  const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

  std::string filename = path + "curved2d.msh";

  std::shared_ptr<Grid> grid = GmshReader<Grid>::read(filename);

  using GridView = typename Grid::LeafGridView;
  GridView gridView = grid->leafGridView();

  // Construct the Brezzi-Douglas-Marini basis
  using Basis = Functions::BrezziDouglasMariniBasis<GridView,1>;
  Basis bdmBasis(gridView);

  // Make a coefficient vector, and fill it with some nontrivial values
  std::vector<double> x(bdmBasis.size());

  for (std::size_t i=0; i<x.size(); i++)
    x[i] = i/M_PI - std::floor(i/M_PI);

  // Generate a Brezzi-Douglas-Marini function with these coefficients
  auto f = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,dim> >(bdmBasis, x);

  passed = passed && Functions::Test::checkNormalContinuity(f);

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e.what() << std::endl;
  return 1;
}
