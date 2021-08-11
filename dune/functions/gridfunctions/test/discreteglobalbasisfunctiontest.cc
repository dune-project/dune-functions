// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/nedelecbasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;

template<class T, int k>
double infinityDiff(const Dune::FieldVector<T,k>& x, const Dune::FieldVector<T,k>& y)
{
  return (x-y).infinity_norm();
}

double infinityDiff(const double& x, const double& y)
{
  return std::fabs(x-y);
}

double infinityDiff(const bool& x, const bool& y)
{
  return std::fabs(x-y);
}

template<class R, class B, class C>
bool checkInterpolationConsistency(B&& basis, C&& x)
{
  using Coeff = typename std::decay<C>::type;
  using Range = R;

  // generate a discrete function
  auto f = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis, x);

  Coeff y;
  interpolate(basis, y, f);
  for (typename std::decay_t<C>::size_type i=0; i<x.size(); ++i)
  {
    if (infinityDiff(x[i],y[i]) > 1e-10)
    {
      std::cout << "Interpolation of DiscreteGlobalBasisFunction differs from original coefficient vector" << std::endl;
      return false;
    }
  }
  return true;
}


int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
  bool passed = true;

  // Generate grid for testing
  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  GridType grid(l,elements);

  using namespace Functions::BasisBuilder;

  const auto& gridView = grid.leafGridView();

  // scalar Lagrange basis with vector valued coefficients
  {
    using Range = FieldVector<double,5>;
    auto f = [](const auto& x){
      Range y;
      for (typename Range::size_type i=0; i<y.size(); ++i)
        y[i] = x[0]+i;
      return y;
    };

    std::vector<Range> x;
    auto feBasis = makeBasis(gridView,lagrange<1>());
    interpolate(feBasis, x, f);
    auto passedThisTest = checkInterpolationConsistency<Range>(feBasis, x);
    std::cout << "checkInterpolationConsistency for scalar Lagrange basis with vector-valued coefficients " << (passedThisTest? " " : " NOT  ") << "successful." << std::endl;
    passed = passed and passedThisTest;
  }

  // scalar Lagrange Basis
  {
    auto feBasis = makeBasis(gridView,lagrange<1>());

    auto f = [](const auto& x){
      return (x.two_norm()<0.5);
    };
    std::vector<bool> x;
    interpolate(feBasis, x, f);
    using Range = bool;
    auto passedThisTest = checkInterpolationConsistency<Range>(feBasis, x);
    std::cout << "checkInterpolationConsistency for scalar Lagrange basis" << (passedThisTest? " " : " NOT ") << "successful." << std::endl;
    passed = passed and passedThisTest;
  }

  // power Lagrange basis
  {
    auto feBasis = makeBasis(gridView,power<2>(lagrange<2>()));
    using Range = FieldVector<double,2>;

    // f(x,y) = (y,x)
    auto f = [](const auto& x){
      return Range{ x[1], x[0] };
    };
    std::vector<Range> x;
    interpolate(feBasis, x, f);
    auto passedThisTest = checkInterpolationConsistency<Range>(feBasis, x);
    std::cout << "checkInterpolationConsistency for power Lagrange basis" << (passedThisTest? " " : " NOT  ") << "successful." << std::endl;
    passed = passed and passedThisTest;
  }

  // Taylor-Hood basis
  {
    auto taylorHoodBasis = makeBasis(
        gridView,
        composite(
          power<dim>(
            lagrange<2>(),
            flatLexicographic()),
          lagrange<1>(),
          flatLexicographic()
          ));
    using namespace Dune::Indices;
    // check with velocity and pressure subspace basis
    {
      auto feBasis = Dune::Functions::subspaceBasis(taylorHoodBasis, _0);

      using Range = FieldVector<double,dim>;

      auto f = [](const auto& x){
        return Range{ x[1], x[0] };
      };
      std::vector<double> x;
      interpolate(feBasis, x, f);
      auto passedThisTest = checkInterpolationConsistency<Range>(feBasis, x);
      std::cout << "checkInterpolationConsistency for velocity part of Taylor-Hood basis" << (passedThisTest? " " : " NOT  ") << "successful." << std::endl;
      passed = passed and passedThisTest;
    }
    {
      auto feBasis = Dune::Functions::subspaceBasis(taylorHoodBasis, _1);

      using Range = double;

      auto f = [](const auto& x){
        return Range{ x[0] };
      };
      std::vector<double> x;
      interpolate(feBasis, x, f);
      auto passedThisTest = checkInterpolationConsistency<Range>(feBasis, x);
      std::cout << "checkInterpolationConsistency for pressure part of Taylor-Hood basis" << (passedThisTest? " " : " NOT  ") << "successful." << std::endl;
      passed = passed and passedThisTest;
    }
  }

  // Raviart-Thomas basis
  {
    auto feBasis = makeBasis(gridView, raviartThomas<0>());
    // coefficients and range of the FE are different here!
    using Coeff = FieldVector<double,1>;
    using Range = FieldVector<double,2>;

    // f(x,y) = (y,x)
    auto f = [](const auto& x){
      return Range{ x[1], x[0] };
    };
    std::vector<Coeff> x;
    interpolate(feBasis, x, f);
    auto passedThisTest = checkInterpolationConsistency<Range>(feBasis, x);
    std::cout << "checkInterpolationConsistency for Raviart-Thomas basis " << (passedThisTest? " " : " NOT  ") << "successful." << std::endl;
    passed = passed and passedThisTest;
  }

  // Nédélec basis
  {
    auto feBasis = makeBasis(gridView, nedelec<1,1>());
    // coefficients and range of the FE are different here!
    using Coeff = FieldVector<double,1>;
    using Range = FieldVector<double,2>;
    // f(x,y) = (y,x)
    auto f = [](const auto& x){
      return Range{ x[1], x[0] };
    };
    std::vector<Coeff> x;
    interpolate(feBasis, x, f);
    auto passedThisTest = checkInterpolationConsistency<Range>(feBasis, x);
    std::cout << "checkInterpolationConsistency for Nédélec basis " << (passedThisTest? " " : " NOT ") << "successful." << std::endl;
    passed = passed and passedThisTest;
  }


  // Sample the function f(x,y) = x on the grid vertices
  // If we use that as the coefficients of a finite element function,
  // we know its integral and can check whether quadrature returns
  // the correct result. Notice that resizing is done by the interpolate method.
  auto feBasis = makeBasis(gridView,lagrange<1>());
  std::vector<FieldVector<double,1> > x;
  auto fAnalytic = [](const auto& x){ return x[0];};
  interpolate(feBasis, x, fAnalytic);
  using Range = FieldVector<double,1>;

  passed = passed and checkInterpolationConsistency<Range>(feBasis, x);

  // generate a discrete function to evaluate the integral
  auto f = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(feBasis, x);

  double exactIntegral = 0.5;

  std::cout << "Testing with raw DiscreteGlobalBasisFunction" << std::endl;
  passed = passed and Dune::Functions::Test::checkGridViewFunction(gridView, f, exactIntegral);

  if (passed)
    std::cout << "All tests passed" << std::endl;

  return passed ? 0: 1;

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
