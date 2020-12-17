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

using namespace Dune;

template <class Range, class Basis>
bool checkBasisConsistency(Basis& basis)
{
  // Define random FE function
  BlockVector<FieldVector<double, 1> > x_rand(basis.size());
  x_rand = 0;
  srand((unsigned int)time(NULL));
  for (auto&& e : x_rand)
    e = double(rand())/double((RAND_MAX));
  using Functions::istlVectorBackend;
  auto globalRandomDiscreteFunction = Functions::makeDiscreteGlobalBasisFunction<Range>(basis, istlVectorBackend(x_rand));

  // Interpolate the random FE function
  BlockVector<FieldVector<double, 1> > x(basis.size());
  interpolate(basis, x, globalRandomDiscreteFunction);

  // Compare both coefficient vectors - they should be the same
  auto diff = x_rand;
  diff -= x;
  return (diff.two_norm() < 1e-12);
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

  auto scalarLagrangeBasis = makeBasis(gridView, lagrange<1>());
  auto passedScalarLagrangeBasisTest = checkBasisConsistency<ScalarRange>(scalarLagrangeBasis);
  std::cout << "Consistency check for scalar Lagrange Basis " << (passedScalarLagrangeBasisTest? "" : "NOT ") << "successful." << std::endl;
  passed = passed and passedScalarLagrangeBasisTest;

  // Test interpolation of constant vector valued function for ...
  using VectorRange = FieldVector<ScalarRange, 2>;

  // ... power Lagrange basis
  auto vectorialLagrangeBasis = makeBasis(gridView, power<2>(lagrange<1>(), flatInterleaved()));
  auto passedVectorialLagrangeBasisTest = checkBasisConsistency<VectorRange>(vectorialLagrangeBasis);
  std::cout << "Consistency check for vectorial Lagrange basis " << (passedVectorialLagrangeBasisTest? "" : "NOT ") << "successful." << std::endl;
  passed = passed and passedVectorialLagrangeBasisTest;

  // ... Raviart-Thomas basis
  auto raviartThomasBasis = makeBasis(gridView, raviartThomas<0>());
  auto passedRaviartThomasBasisTest = checkBasisConsistency<VectorRange>(raviartThomasBasis);
  std::cout << "Consistency check for Raviart-Thomas basis " << (passedRaviartThomasBasisTest? "" : "NOT ") << "successful." << std::endl;
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
