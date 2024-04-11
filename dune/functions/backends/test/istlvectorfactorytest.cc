// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/backends/istlvectorfactory.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>


namespace CD = Dune::Functions::ContainerDescriptors;

void checkVectorFactory (Dune::TestSuite& test)
{
  using namespace Dune::Indices;

  auto vec0 = Dune::Functions::istlVectorFactory<double>(CD::Value{});
  static_assert(std::is_same_v<decltype(vec0),double>);

  auto vec1 = Dune::Functions::istlVectorFactory<double>(CD::FlatArray<2>{});
  static_assert(std::is_same_v<decltype(vec1),Dune::FieldVector<double,2>>);

  auto vec2 = Dune::Functions::istlVectorFactory<double>(CD::FlatVector{10});
  static_assert(std::is_same_v<decltype(vec2),Dune::BlockVector<double>>);
  test.check(vec2.size() == 10, "vec2.size");

  auto vec3 = Dune::Functions::istlVectorFactory<double>(CD::UniformVector<CD::FlatArray<2>>{10});
  static_assert(std::is_same_v<decltype(vec3),Dune::BlockVector<Dune::FieldVector<double,2>>>);
  test.check(vec3.size() == 10, "vec3.size");

  // more complicated test
  CD::Tuple<CD::Array<CD::FlatVector,3>,CD::FlatVector> stokes{
    CD::Array<CD::FlatVector,3>{CD::FlatVector{10},CD::FlatVector{10},CD::FlatVector{10}}, CD::FlatVector{5} };

  auto vec4 = Dune::Functions::istlVectorFactory<double>(stokes);
  test.check(vec4.size() == 2, "vec4.size");
  test.check(vec4[_0].size() == 3, "vec4[0].size");
  test.check(vec4[_1].size() == 5, "vec4[01.size");
  test.check(vec4[_0][0].size() == 10, "vec4[0][0].size");
  test.check(vec4[_0][1].size() == 10, "vec4[0][1].size");
  test.check(vec4[_0][2].size() == 10, "vec4[0][2].size");

  auto vec5 = Dune::Functions::istlVectorFactory<double>(stokes[_0]);
  test.check(vec5.size() == vec4[_0].size(), "vec5 == vec4[0]");
}

int main (int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;
  checkVectorFactory(test);

  return test.exit();
}
