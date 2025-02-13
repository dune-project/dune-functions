// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#include <config.h>

#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/backends/containerfactory.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>


namespace CD = Dune::Functions::ContainerDescriptors;

void checkContainerFactory (Dune::TestSuite& test)
{
  using namespace Dune::Indices;

  auto vec0 = Dune::Functions::makeContainer(CD::Value{}, 42);
  static_assert(std::is_same_v<decltype(vec0),int>);
  test.check(vec0==42, "vec0.value");

  auto vec1 = Dune::Functions::makeContainer<int>(CD::FlatArray<2>{});
  static_assert(std::is_same_v<decltype(vec1),std::array<int,2>>);

  auto vec2 = Dune::Functions::makeContainer<int>(CD::FlatVector{10});
  static_assert(std::is_same_v<decltype(vec2),std::vector<int>>);
  test.check(vec2.size() == 10, "vec2.size");

  auto vec3 = Dune::Functions::makeContainer<int>(CD::UniformVector<CD::FlatArray<2>>{10});
  static_assert(std::is_same_v<decltype(vec3),std::vector<std::array<int,2>>>);
  test.check(vec3.size() == 10, "vec3.size");

  // more complicated test
  CD::Tuple<CD::Array<CD::FlatVector,3>,CD::FlatVector> stokes{
    CD::Array<CD::FlatVector,3>{CD::FlatVector{10},CD::FlatVector{10},CD::FlatVector{10}}, CD::FlatVector{5} };

  auto vec4 = Dune::Functions::makeContainer<int>(stokes);
  static_assert(std::is_same_v<decltype(vec4), Dune::TupleVector<std::array<std::vector<int>, 3>, std::vector<int>>>);
  test.check(vec4.size() == 2, "vec4.size");
  test.check(vec4[_0].size() == 3, "vec4[0].size");
  test.check(vec4[_1].size() == 5, "vec4[01.size");
  test.check(vec4[_0][0].size() == 10, "vec4[0][0].size");
  test.check(vec4[_0][1].size() == 10, "vec4[0][1].size");
  test.check(vec4[_0][2].size() == 10, "vec4[0][2].size");

  auto vec5 = Dune::Functions::makeContainer<int>(stokes[_0]);
  test.check(vec5.size() == vec4[_0].size(), "vec5 == vec4[0]");
  static_assert(std::is_same_v<decltype(vec5), std::array<std::vector<int>, 3>>);
  test.check(vec4.size() == 2, "vec4.size");
}

int main (int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;
  checkContainerFactory(test);

  return test.exit();
}
