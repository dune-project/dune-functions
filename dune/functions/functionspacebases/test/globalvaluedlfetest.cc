// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define  DUNE_FUNCTION_HH_SILENCE_DEPRECATION

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/functions/functionspacebases/nedelecbasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>


using namespace Dune;


template<class Basis>
void checkBasisFEs(const Basis& basis) {
  auto localView = basis.localView();
  for (const auto& element : elements(basis.gridView()))
  {
    localView.bind(element);
    testFE(localView.tree().finiteElement(), DisableNone, 1 /* diffOrder */);
  }
}


int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  using namespace Functions::BasisFactory;

  // Check with UGGrid
  {
    // Generate grid for testing
    const int dim = 2;
    using Grid = UGGrid<dim>;

    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    std::string filename = path + "curved2d.msh";
    std::shared_ptr<Grid> grid = GmshReader<Grid>::read(filename);

    auto gridView = grid->leafGridView();

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(div)-conforming space
    //  We use the Raviart-Thomas basis.
    ///////////////////////////////////////////////////////////////////////

    checkBasisFEs(makeBasis(gridView, raviartThomas<0>()));
    checkBasisFEs(makeBasis(gridView, raviartThomas<1>()));

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(curl)-conforming space
    //  We use the Nedelec basis of the first kind.
    ///////////////////////////////////////////////////////////////////////

    checkBasisFEs(makeBasis(gridView, nedelec<1,1,double>()));
  }

  // Check with YaspGrid
  {
    // Generate grid for testing
    auto grid = YaspGrid<2>({1.0, 1.0}, {5,5});
    auto gridView = grid.leafGridView();

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(div)-conforming space
    //  We use the Raviart-Thomas basis.
    ///////////////////////////////////////////////////////////////////////

    checkBasisFEs(makeBasis(gridView, raviartThomas<0>()));
    checkBasisFEs(makeBasis(gridView, raviartThomas<1>()));
    checkBasisFEs(makeBasis(gridView, raviartThomas<2>()));

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(curl)-conforming space
    //  We use the Nedelec basis of the first kind.
    ///////////////////////////////////////////////////////////////////////

    checkBasisFEs(makeBasis(gridView, nedelec<1,1,double>()));
  }

} catch (Exception &e)
{
  std::cerr << "Dune reported error: " << e.what() << std::endl;
  return 1;
}
