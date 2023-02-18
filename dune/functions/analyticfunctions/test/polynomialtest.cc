// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/functions/analyticfunctions/polynomial.hh>



// Compare coefficients of different polynomials
template<class P0, class P1>
bool polynomialsEqual(const P0& p0, const P1& p1) {
  std::vector<double> c0;
  Dune::Hybrid::forEach(p0.coefficients(), [&](auto z) {
    c0.push_back(z);
  });
  std::vector<double> c1;
  Dune::Hybrid::forEach(p1.coefficients(), [&](auto z) {
    c1.push_back(z);
  });
  auto n0 = c0.size();
  auto n1 = c1.size();
  for([[maybe_unused]] auto i : Dune::range(n0,std::max(n0, n1)))
    c0.push_back(0);
  for([[maybe_unused]] auto i : Dune::range(n1,std::max(n0, n1)))
    c1.push_back(0);
  return (c0==c1);
}



using Dune::Functions::Polynomial;
using Dune::Functions::makePolynomial;



template<class... Args>
Dune::TestSuite testDynamicPolynomial()
{
  Dune::TestSuite suite;

  // Check construction from initializer list
  using P = Polynomial<Args...>;
  auto p = P({ 1, 2, 3, -4});

  {
    double x = 42;
    double y = 1.0 + 2.0*x + 3.0*x*x - 4.0*x*x*x;
    suite.check(std::fabs(p(x) - y) < 1e-14);
  }

  // Check construction from l-value reference
  using Coeff = typename P::Coefficients;
  auto c = Coeff({1, 2, 3, -4});
  suite.check(polynomialsEqual(p, P(c)));

  // Check default construction
  {
    [[maybe_unused]] P p1;
    [[maybe_unused]] auto p2 = P();
  }

  // Check copy construction
  {
    auto p2 = p;
    suite.check(polynomialsEqual(p, p2));
    P p3(p2);
    suite.check(polynomialsEqual(p, p3));
  }

  // Check copy assignment
  {
    auto p1 = P({2,3,4});
    auto p2 = P({3,4,5});
    suite.check(not polynomialsEqual(p1, p2));
    p1 = p2;
    suite.check(polynomialsEqual(p1, p2));
  }

  // Check move assignment
  {
    auto p1 = P({2,3,4});
    auto p2 = P({3,4,5});
    p1 = P({3,4,5});
    suite.check(polynomialsEqual(p1, p2));
  }

  // Check some value
  suite.check(p(0) == 1);

  // Check derivatives
  {
    auto dp = P({ 2, 6, -12 });
    auto ddp = P({ 6, -24 });
    auto dddp = P({ -24 });
    auto zero = P();

    suite.check(polynomialsEqual(derivative(p), dp));
    suite.check(polynomialsEqual(derivative(derivative(p)), ddp));
    suite.check(polynomialsEqual(derivative(derivative(derivative(p))), dddp));
    suite.check(polynomialsEqual(derivative(derivative(derivative(derivative(p)))), zero));
    suite.check(polynomialsEqual(derivative(derivative(derivative(derivative(derivative(p))))), zero));
  }

  return suite;
}



template<class K>
Dune::TestSuite testIntegerSequencePolynomial()
{
  Dune::TestSuite suite;

  auto p = Polynomial<K, std::integer_sequence<int, 1,2,3,-4>>();
  auto dp = Polynomial<K, std::integer_sequence<int,2,6,-12>>();
  auto ddp = Polynomial<K, std::integer_sequence<int,6, -24>>();
  auto dddp = Polynomial<K, std::integer_sequence<int,-24>>();
  auto zero = Polynomial<K, std::integer_sequence<int>>();

  suite.check(derivative(p)==dp);
  suite.check(derivative(derivative(p))==ddp);
  suite.check(derivative(derivative(derivative(p)))==dddp);
  suite.check(derivative(derivative(derivative(derivative(p))))==zero);
  suite.check(derivative(derivative(derivative(derivative(derivative(p)))))==zero);

  suite.check(polynomialsEqual(p, Polynomial<K>({1,2,3,-4})));
  suite.check(polynomialsEqual(derivative(p), Polynomial<K>({2,6,-12})));
  suite.check(polynomialsEqual(derivative(derivative(p)), Polynomial<K>({6,-24})));
  suite.check(polynomialsEqual(derivative(derivative(derivative(p))), Polynomial<K>({-24})));
  suite.check(polynomialsEqual(derivative(derivative(derivative(derivative(p)))), Polynomial<K>()));
  suite.check(polynomialsEqual(derivative(derivative(derivative(derivative(derivative(p))))), Polynomial<K>()));

  {
    double x = 42;
    double y = 1.0 + 2.0*x + 3.0*x*x - 4.0*x*x*x;
    suite.check(std::fabs(p(x) - y) < 1e-14);
  }

  return suite;
}



template<class K>
Dune::TestSuite testTuplePolynomial()
{
  Dune::TestSuite suite;

  using namespace Dune::Indices;

  auto p = makePolynomial<K>(std::tuple(1ul, 2ul, _3, -4l));
  auto dp = makePolynomial<K>(std::tuple(2ul, _6, -12l));
  auto ddp = makePolynomial<K>(std::tuple(_6, -24l));
  auto dddp = makePolynomial<K>(std::tuple(-24l));
  auto zero = makePolynomial<K>(std::tuple());

  suite.check(derivative(p)==dp);
  suite.check(derivative(derivative(p))==ddp);
  suite.check(derivative(derivative(derivative(p)))==dddp);
  suite.check(derivative(derivative(derivative(derivative(p))))==zero);
  suite.check(derivative(derivative(derivative(derivative(derivative(p)))))==zero);

  suite.check(polynomialsEqual(p, Polynomial<K>({1,2,3,-4})));
  suite.check(polynomialsEqual(derivative(p), Polynomial<K>({2,6,-12})));
  suite.check(polynomialsEqual(derivative(derivative(p)), Polynomial<K>({6,-24})));
  suite.check(polynomialsEqual(derivative(derivative(derivative(p))), Polynomial<K>({-24})));
  suite.check(polynomialsEqual(derivative(derivative(derivative(derivative(p)))), Polynomial<K>()));
  suite.check(polynomialsEqual(derivative(derivative(derivative(derivative(derivative(p))))), Polynomial<K>()));

  {
    double x = 42;
    double y = 1.0 + 2.0*x + 3.0*x*x -4.0*x*x*x;
    suite.check(std::fabs(p(x) - y) < 1e-14);
  }

  return suite;
}



int main(int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite suite;

  suite.subTest(testDynamicPolynomial<int>());
  suite.subTest(testDynamicPolynomial<double>());
  suite.subTest(testDynamicPolynomial<double, std::array<double,4>>());
  suite.subTest(testDynamicPolynomial<double, std::vector<double>>());

  suite.subTest(testIntegerSequencePolynomial<int>());
  suite.subTest(testIntegerSequencePolynomial<double>());

  suite.subTest(testTuplePolynomial<int>());
  suite.subTest(testTuplePolynomial<double>());

  return suite.exit();
}
