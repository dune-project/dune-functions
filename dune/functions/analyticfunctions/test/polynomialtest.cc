// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/functions/analyticfunctions/polynomial.hh>

using namespace Dune;
using namespace Dune::Functions;

TestSuite
testPolynomial()
{
  TestSuite suite;

  auto polynomialsEqual = [](const auto& p_0, const auto& p_1) {
    return p_0.coefficients() == p_1.coefficients();
  };

  using P = Polynomial<int>;
  auto p = P({ 1, 2, 3 });

  suite.check(p(0) == 1);

  // check different copy construction variants to see if they compile:
  auto p2 = p;
  P p3(p2);

  // check if they have the same coefficients
  suite.check(polynomialsEqual(p, p2));
  suite.check(polynomialsEqual(p, p3));

  // Use assignment operators
  p2 = p3;
  suite.check(polynomialsEqual(p2, p3));

  p3 = std::move(p2);
  suite.check(polynomialsEqual(p, p3));

  return suite;
}

int
main(int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  TestSuite suite;
  suite.subTest(testPolynomial());

  return suite.exit();
}
