// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define  DUNE_FUNCTION_HH_SILENCE_DEPRECATION

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/functions/functionspacebases/raviartthomasbasis.hh>


using namespace Dune;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  // Generate grid for testing
  const int dim = 2;
  using Grid = UGGrid<dim>;

  const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

  std::string filename = path + "curved2d.msh";

#if 0
  std::shared_ptr<Grid> grid = GmshReader<Grid>::read(filename);
#else
  GridFactory<Grid> factory;

  factory.insertVertex({0,0});
  factory.insertVertex({0.5,0});
  factory.insertVertex({0.5 ,0.5});

  factory.insertElement(GeometryTypes::simplex(dim), {0,1,2});
  std::shared_ptr<Grid> grid = factory.createGrid();
#endif

  using GridView = typename Grid::LeafGridView;
  GridView gridView = grid->leafGridView();

  ///////////////////////////////////////////////////////////////////////
  //  Test GlobalValuedLocalFiniteElement for a H(div)-conforming space
  //  We use the Raviart-Thomas basis.
  ///////////////////////////////////////////////////////////////////////

  using Basis = Functions::RaviartThomasBasis<GridView,0>;
  Basis rtBasis(gridView);

  auto localView = rtBasis.localView();

  for (const auto& element : elements(gridView))
  {
    localView.bind(element);

    testFE(localView.tree().finiteElement(), DisableNone, 1 /* diffOrder */);
  }

} catch (Exception &e)
{
  std::cerr << "Dune reported error: " << e.what() << std::endl;
  return 1;
}
