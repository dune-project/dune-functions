// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>

#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/bvector.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
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
  std::array<int,dim> N = {1,1};
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

  // { definition_basis1_begin }
  Functions::PQkNodalBasis<GridView,2> basis1(gridView);
  // { definition_basis1_end }

  // { definition_x1_begin }
  std::vector<double> x1;
  // { definition_x1_end }

  // { interpolation1_begin }
  interpolate(basis1, x1, f1);
  // { interpolation1_end }

  // Interpolate a vector-valued function using a vector-valued basis
  // { interpolation2_begin }
  auto f2 = [](const FieldVector<double,2>& x)
  {
    return x;
  };

  //Functions::PQkNodalBasis<GridView,2> basis2(gridView);
  using namespace Functions::BasisBuilder;

  auto basis2 = makeBasis(
    gridView,
      power<dim>(
        pq<2>()
      )
  );

  BlockVector<FieldVector<double,2> > x2;

  interpolate(basis2, x2, f2);
  // { interpolation2_end }

  // Interpolate a vector-valued function using a scalar basis
  // { interpolation3_begin }
  auto f3 = [](const FieldVector<double,dim>& x) { return x; };

  Functions::PQkNodalBasis<GridView,2> basis3(gridView);
  BlockVector<FieldVector<double,3> > x3;
  interpolate(basis3, x3, f3);
  // { interpolation3_end }

  // { interpolation4_begin }
  // **********************************************
  //  FEHLT NOCH!
  // **********************************************
  // { interpolation4_end }

  // { interpolation5_begin }
  // **********************************************
  //  FEHLT NOCH!
  // **********************************************
#if 0
  BitSetVector<3> boundary = ...;  // initialize somehow

  for (size_t i=0; i<boundary.size(); i++)
    if ( /* i is boundary */ )
      boundary[i] = "011";   // works like std::bitset

  interpolate...
#endif
  // { interpolation5_end }

  return 0;
}
