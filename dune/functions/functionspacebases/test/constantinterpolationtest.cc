// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/bvector.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

using namespace Dune;

template <class Range, class Basis, class InputFunction, class GridView>
bool interpolateAndCompareFunctions(Basis& basis, InputFunction& inputFunction, GridView& gridView)
{
  // Interpolate
  BlockVector<FieldVector<double, 1> > x(basis.size());
  x = 0;
  interpolate(basis, x, inputFunction);
  using Functions::istlVectorBackend;
  auto globalDiscreteFunction = Functions::makeDiscreteGlobalBasisFunction<Range>(basis, istlVectorBackend(x));
  auto localDiscreteFunction = localFunction(globalDiscreteFunction);
  auto localInputFunction = localFunction(Functions::makeGridViewFunction(inputFunction, gridView));

  // Compare - Since a constant function is used, evaluate and compare functions only in a single point.
  bool functionsEqual = true;
  for (const auto& element: elements(gridView))
  {
    localInputFunction.bind(element);
    localDiscreteFunction.bind(element);
    const auto geometry = element.geometry();
    const auto center = geometry.local(geometry.center());
    auto diff = localInputFunction(center);
    diff -= localDiscreteFunction(center);
    functionsEqual = functionsEqual and (std::sqrt(diff * diff) < 1e-6);
  }

  return functionsEqual;
}

int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);
  bool passed = true;

  // Generate grid for testing
  const int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> l(1.0);
  std::array<int,dim> elements = {{2, 2}};
  auto gridPtr = std::make_unique<Grid>(l, elements);
  auto& grid = *gridPtr;
  grid.globalRefine(2);
  auto gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////
  using namespace Functions::BasisBuilder;

  // Test interpolation of constant function for
  // ... scalar Lagrange basis
  using ScalarRange = double;
  auto scalarValuedFunction = [](const auto& x) { return 2;};

  auto scalarLagrangeBasis = makeBasis(gridView, lagrange<1>());
  auto passedScalarLagrangeBasisTest = interpolateAndCompareFunctions<ScalarRange>(scalarLagrangeBasis, scalarValuedFunction, gridView);
  std::cout << "Constant interpolation for scalar Lagrange Basis " << (passedScalarLagrangeBasisTest? "" : "NOT ") << "successful." << std::endl;
  passed = passed and passedScalarLagrangeBasisTest;

  // Test interpolation of constant vector valued function for ...
  using VectorRange = FieldVector<ScalarRange, 2>;
  auto vectorValuedFunction = [](const auto& x) { return FieldVector<double,2>{1, 2};};

  // ... power Lagrange basis
  auto vectorialLagrangeBasis = makeBasis(gridView, power<2>(lagrange<1>(), flatInterleaved()));
  auto passedVectorialLagrangeBasisTest = interpolateAndCompareFunctions<VectorRange>(vectorialLagrangeBasis, vectorValuedFunction, gridView);
  std::cout << "Constant interpolation for vectorial Lagrange basis " << (passedVectorialLagrangeBasisTest? "" : "NOT ") << "successful." << std::endl;
  passed = passed and passedVectorialLagrangeBasisTest;

  // ... Raviart-Thomas basis
  auto raviartThomasBasis = makeBasis(gridView, raviartThomas<0>());
  auto passedRaviartThomasBasisTest = interpolateAndCompareFunctions<VectorRange>(raviartThomasBasis, vectorValuedFunction, gridView);
  std::cout << "Constant interpolation for Raviart-Thomas basis " << (passedRaviartThomasBasisTest? "" : "NOT ") << "successful." << std::endl;
  passed = passed and passedRaviartThomasBasisTest;

  if (passed)
    std::cout << "All tests passed" << std::endl;

  return passed ? 0 : 1;

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
