// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <array>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/functions/analyticfunctions/monomialset.hh>




// Create n^dim uniformly distributed points in [0,1]^dim
template<class F, unsigned int dim>
std::vector<Dune::FieldVector<F, dim>> samplePoints(unsigned int n)
{
  unsigned long size = std::pow(n, dim);
  auto points = std::vector<Dune::FieldVector<F, dim>>(size);
  if constexpr (dim==1)
  {
    for(auto k : Dune::range(n))
      points[k] = F(k)/F(n-1);
  }
  else
  {
    auto points_low = samplePoints<F, dim-1>(n);
    auto n_low = points_low.size();
    for(auto k : Dune::range(size))
    {
      auto k_0 = k % n_low;
      auto k_1 = k / n_low;
      for(auto i : Dune::range(dim-1))
        points[k][i] = points_low[k_0][i];
      points[k][dim-1] = F(k_1)/F(n-1);
    }
  }
  return points;
}



template<class F, unsigned long dim, unsigned long maxOrder>
Dune::TestSuite testMonomialSet(unsigned long n, double tol)
{
  Dune::TestSuite suite;

  auto monomials = Dune::Functions::MonomialSet<F, dim, maxOrder>{};

  auto p = monomials.exponents();
  auto size = p.size();

  for(auto x : samplePoints<F, dim>(n))
  {
    auto y = monomials(x);
    suite.check(y.size() == size);
    for(auto i : Dune::range(size))
    {
      auto yy = F(1.0);
      for(auto k : Dune::range(dim))
        yy *= std::pow(x[k], p[i][k]);
      suite.check(std::fabs(y[i] - yy) < tol)
        << "Monomial(dim=" << dim << ",maxOrder=" << maxOrder << ",index=" << i << ") value incorrect";
    }
  }

  auto D_monomials = derivative(monomials);
  for(auto x : samplePoints<F, dim>(n))
  {
    auto y = D_monomials(x);
    suite.check(y.N() == size);
    suite.check(y.M() == dim);
    for(auto i : Dune::range(size))
    {
      for(auto j : Dune::range(dim))
      {
        auto yy = F(1.0);
        for(auto k : Dune::range(dim))
        {
          if (p[i][k]-(k==j) > 0)
            yy *= std::pow(x[k], p[i][k]-(k==j));
        }
        yy *= p[i][j];
        suite.check(std::fabs(y[i][j] - yy) < tol)
          << "Monomial(dim=" << dim << ",maxOrder=" << maxOrder << ",index=" << i << ") derivative incorrect";
      }
    }
  }

  auto H_monomials = derivative(derivative(monomials));
  for(auto x : samplePoints<F, dim>(n))
  {
    auto y = H_monomials(x);
    suite.check(y.size() == size)<<"Wrong Size";
    suite.check(y[0].M() == dim)<<"Wrong M";
    suite.check(y[0].N() == dim)<<"Wrong N";

    for(auto i : Dune::range(size))
    {
      for(auto j : Dune::range(dim))
      {
        for (auto l: Dune::range(dim)){
          auto yy = F(1.0);
          for(auto k : Dune::range(dim))
          {
            if (p[i][k]-(k==j)-(k==l) > 0)
              yy *= std::pow(x[k], p[i][k] - int(k == j) - int(k == l));
          }
          if (j == l)
            yy *= p[i][j] * (int(p[i][j]) - 1.);
          else
            yy *= p[i][j] * p[i][l];

          suite.check(std::fabs(y[i][j][l] - yy) < tol)
            << "Monomial(dim=" << dim << ",maxOrder=" << maxOrder << ",index=" << i << ", exponents= "<< p[i][0] << p[i][1]<<") hessian component ["<<l<<","<<j<<"] was "<<y[i][j][l]<<" but "<<yy<<" was expected.";
        }
      }
    }
  }

  return suite;
}



int main(int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite suite;

  suite.subTest(testMonomialSet<double, 1, 0>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 1, 1>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 1, 2>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 1, 3>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 1, 4>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 1, 5>(10, 1e-14));

  suite.subTest(testMonomialSet<double, 2, 0>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 2, 1>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 2, 2>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 2, 3>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 2, 4>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 2, 5>(10, 1e-14));

  suite.subTest(testMonomialSet<double, 3, 0>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 3, 1>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 3, 2>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 3, 3>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 3, 4>(10, 1e-14));
  suite.subTest(testMonomialSet<double, 3, 5>(10, 1e-14));
  return suite.exit();
}
