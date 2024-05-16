// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_INTERPOLATETEST_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_INTERPOLATETEST_HH

#include <tuple>
#include <utility>

#include <dune/common/test/testsuite.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/functionspacebases/interpolate.hh>


double infinityDiff(const double& x, const double& y)
{
  return std::fabs(x-y);
}

double infinityDiff(const bool& x, const bool& y)
{
  return std::fabs(x-y);
}

template<class X, class Y>
double infinityDiff(const X& x, const Y& y)
{
  if (x.size() != y.size())
    return false;
  double diff = 0;
  Dune::Hybrid::forEach(Dune::range(Dune::Hybrid::size(x)), [&](auto i) {
    auto&& xi = Dune::Hybrid::elementAt(x, i);
    auto&& yi = Dune::Hybrid::elementAt(y, i);
    diff = std::max(diff, infinityDiff(xi, yi));
  });
  return diff;
}

template<class Range, class Basis, class C>
Dune::TestSuite checkInterpolateConsistency(Basis basis, C&& x)
{
  using Coefficients = std::decay_t<C>;

  Dune::TestSuite suite("interpolate consistency check");
  double coeffTol = 1e-10;

  // generate a discrete function
  auto fGridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis, x);

  // Check both functions
  {
    const auto& f = fGridFunction;
    Coefficients y;
    Dune::Functions::interpolate(basis, y, f);

    suite.check(infinityDiff(x, y) < coeffTol)
      << "Interpolation of DiscreteGlobalBasisFunction via local operator() differs from original coefficient vector" << std::endl;
  }

  return suite;
}



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_INTERPOLATETEST_HH
