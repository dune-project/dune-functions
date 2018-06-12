// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/tuplevector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

using namespace Dune;


int main (int argc, char *argv[])
{
  //Maybe initialize Mpi
  MPIHelper::instance(argc, argv);

  // make grid
  const int dim = 2;
  FieldVector<double,dim> L(1.0);
  std::array<int,dim> N = {{1,1}};
  std::bitset<dim> periodic(false);
  YaspGrid<2> grid(L,N);
  using GridView = YaspGrid<2>::LeafGridView;
  GridView gridView = grid.leafGridView();

  // Create closed-form function to interpolate
  // { definition_f1_begin }
  auto f1 = [](const FieldVector<double,2>& x)
  {
    return exp(-1.0*x.two_norm2());
  };
  // { definition_f1_end }

  // { definition_p2basis_begin }
  Functions::LagrangeBasis<GridView,2> p2basis(gridView);
  // { definition_p2basis_end }

  // { definition_x1_begin }
  std::vector<double> x1;
  // { definition_x1_end }

  // { interpolation1_begin }
  interpolate(p2basis, x1, f1);
  // { interpolation1_end }

  // Interpolate a vector-valued function using a vector-valued basis
  // { taylorhood_basis_begin }
  using namespace Functions::BasisBuilder;

  auto taylorHoodBasis = makeBasis(
    gridView,
    composite(
      power<dim>(
        lagrange<2>(),
        flatLexicographic()),
      lagrange<1>(),
      flatLexicographic()
    ));
  // { taylorhood_basis_end }

  // { taylorhood_vector_begin }
  BlockVector<FieldVector<double,1>> x2;
  // { taylorhood_vector_end }


  // { taylorhood_pressure_begin }
  using namespace Indices;
  interpolate(subspaceBasis(taylorHoodBasis, _1), x2, f1);
  // { taylorhood_pressure_end }

  // { taylorhood_velocity_begin }
  auto f2 = [](const FieldVector<double,2>& x) {
    return x;
  };
  interpolate(subspaceBasis(taylorHoodBasis, _0), x2, f2);
  // { taylorhood_velocity_end }

  // { setup_mask_begin }
  BlockVector<FieldVector<char,1>> isBoundary;
  auto isBoundaryBackend = Functions::istlVectorBackend(isBoundary);
  isBoundaryBackend.resize(taylorHoodBasis);
  isBoundary = false;
  forEachBoundaryDOF(subspaceBasis(taylorHoodBasis, _0),
    [&] (auto&& index) {
      isBoundaryBackend[index] = true;
    });
  // { setup_mask_end }

  // { masked_interpolation_begin }
  interpolate(subspaceBasis(taylorHoodBasis, _1), x2, f2, isBoundary);
  // { masked_interpolation_end }


  return 0;
}
