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

  ///////////////////////////////////
  //   Generate a simple grid
  ///////////////////////////////////
  const int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> l(1.0);
  std::array<int,dim> elements = {{2, 2}};
  auto gridPtr = std::make_unique<Grid>(l, elements);

  auto& grid = *gridPtr;
  grid.globalRefine(2);

  auto gridView = grid.leafGridView();
  using GridView = decltype(gridView);

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////
  using namespace Functions::BasisBuilder;
  using namespace TypeTree::Indices;

  auto vectorialLagrangeBasis = makeBasis(gridView, power<2>(lagrange<1>(), flatInterleaved()));
  auto raviartThomasBasis = makeBasis(gridView, raviartThomas<0>());

  // Test interpolation of constant function for scalar Lagrange basis
  using ScalarRange = double;
  auto scalarValuedFunction = [](const auto& x) { return 2;};

  auto scalarLagrangeBasis = makeBasis(gridView, lagrange<1>());
  auto scalarLagrangeBasisTest = interpolateAndCompareFunctions<ScalarRange>(scalarLagrangeBasis, scalarValuedFunction, gridView);
  std::cout << "Constant interpolation for scalar Lagrange Basis " << (scalarLagrangeBasisTest? "" : "NOT ") << "successful." << std::endl;
  assert(scalarLagrangeBasisTest == true);

  // Test interpolation of constant function for vectorial bases
  using VectorRange = FieldVector<ScalarRange, 2>;
  auto vectorValuedFunction = [](const auto& x) { return FieldVector<double,2>{1, 2};};

  // ... Lagrange basis
  auto vectorialLagrangeBasisTest = interpolateAndCompareFunctions<VectorRange>(vectorialLagrangeBasis, vectorValuedFunction, gridView);
  std::cout << "Constant interpolation for vectorial Lagrange basis " << (vectorialLagrangeBasisTest? "" : "NOT ") << "successful." << std::endl;
  assert(raviartThomasTest == true);

  // ... Raviart-Thomas basis
  auto raviartThomasBasisTest = interpolateAndCompareFunctions<VectorRange>(raviartThomasBasis, vectorValuedFunction, gridView);
  std::cout << "Constant interpolation for Raviart-Thomas basis " << (raviartThomasBasisTest? "" : "NOT ") << "successful." << std::endl;
  assert(raviartThomasTest == true);

 }
// Error handling
 catch (Exception& e) {
    std::cout << e.what() << std::endl;
 }
