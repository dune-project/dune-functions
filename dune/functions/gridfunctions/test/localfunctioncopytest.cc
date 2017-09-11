// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/test/testsuite.hh> // We use exceptions

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>

#include <dune/istl/bvector.hh>
using namespace Dune;

template<class F, class GV>
auto checkCopyConstruction(F&& f, const GV& gv) {
  TestSuite suite;
  auto localF = localFunction(f);
  suite.check(not localF.bound(), "Check if LocalFunction was not marked bound too early.");

  for (const auto& element: elements(gv)) {
    // bind localF
    localF.bind(element);
    suite.check(localF.bound(), "Check if LocalFunction was bound.");

    auto localF_copy = localF;
    suite.check(localF_copy.bound(), "Check if copied LocalFunction was bound.");

    // Test evaluation
    using Domain = typename std::decay_t<F>::Domain;
    auto xval = Domain(0);
    // evaluate the copied localFunction. This gave a segmentation fault earlier
    auto copy_eval = localF_copy(xval);
    auto eval = localF(xval);
    suite.check((copy_eval == eval), "Check if copied and original local function give same result");
  }
  return suite;
}

int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);

  TestSuite suite;

  // Grid
  using GridType =Dune::YaspGrid<2>;
  auto grid = Dune::StructuredGridFactory<GridType>
    ::createCubeGrid({{0,0}}, {{1,1}}, {{2,2}});

  // Basis
  using Basis = Dune::Functions::PQ1NodalBasis<GridType::LeafGridView>;
  Basis basis(grid->leafGridView());

  // some artificial coefficients
  using Vector = Dune::BlockVector<Dune::FieldVector<double, 1> >;
  Vector x;
  x.resize(basis.dimension());
  auto xbe = Dune::Functions::HierarchicVectorWrapper<Vector>(x);
  x=1.0;

  // create a discrete global basis function
  auto f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
    (basis, Dune::TypeTree::hybridTreePath(), xbe);

  // test with DiscreteGlobalBasisFunction
  suite.subTest(checkCopyConstruction(f, grid->leafGridView()));
  return suite.exit();
}
