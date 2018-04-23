// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <array>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/persistentgridview.hh>

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  using namespace Dune;
  using namespace Dune::Functions;
  using namespace Dune::Functions::BasisFactory;

  TestSuite test;

  // Generate grid for testing
  const int dim = 2;
  using Grid = YaspGrid<dim>;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};

  Grid grid(l,elements);


  // Helper to change thrown exception to checkable condition
  auto notThrown = [](auto&& f) {
    try {
      f();
      return true;
    }
    catch (Exception e) {
      return false;
    }
  };

  {
    auto gridView = grid.leafGridView();

    auto persistentBasis = makeBasis(Experimental::persistentGridView(gridView), lagrange<3>());

    // Check if filling the PersistentGridView lazily
    // in mutable state works
    test.require(notThrown([&]() {
      auto localView = persistentBasis.localView();
      for(const auto& e : Dune::elements(gridView))
        localView.bind(e);
      }))
      << "Filling PersistentGridView in mutable state failed.";

    persistentBasis.gridView().indexSet().setMutable(false);

    // Check if accessing the PersistentGridView in non-mutable state works
    test.require(notThrown([&]() {
      auto localView = persistentBasis.localView();
      for(const auto& e : Dune::elements(gridView))
        localView.bind(e);
      }))
      << "Accessing PersistentGridView indices in non-mutable state failed.";

    grid.globalRefine(1);

    // Check if accessing the PersistentGridView in non-mutable state works
    // after grid refinement
    test.require(notThrown([&]() {
      // Create coarse grid view because you cannot iterate over
      // a PersistentGridView
      auto coarseGridView = grid.levelGridView(grid.maxLevel()-1);
      auto localView = persistentBasis.localView();
      for(const auto& e : Dune::elements(coarseGridView))
        localView.bind(e);
      }))
      << "Accessing PersistentGridView indices in non-mutable state failed after grid refinement.";

    // Check if the coarse elements are known to the PersistentGridView
    {
      bool containsCoarseElements = true;
      auto pgv = persistentBasis.gridView();

      auto coarseGridView = grid.levelGridView(grid.maxLevel()-1);
      for (const auto& e : Dune::elements(coarseGridView)) {
        containsCoarseElements = containsCoarseElements && pgv.contains(e);
      }
      test.require(containsCoarseElements, "Check contains()") << "Some elements that are part of the PersistentGridView are not contained.";
    }

    // Reversely, check if the new elements are unknown to the PersistentGridView
    {
      bool containsNoFineElements = true;
      auto pgv = persistentBasis.gridView();

      for (const auto& e : Dune::elements(grid.leafGridView())) {
        containsNoFineElements = containsNoFineElements && not pgv.contains(e);
      }
      test.require(containsNoFineElements, "Check contains()") << "Some elements that are part of the PersistentGridView are not contained.";
    }
  }

  return test.exit();
}
