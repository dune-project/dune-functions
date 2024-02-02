// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <type_traits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>


namespace CD = Dune::Functions::ContainerDescriptors;


void checkHierarchic (Dune::TestSuite& test)
{
  using namespace Dune::Indices;
  CD::Tuple<CD::Array<CD::FlatVector,3>,CD::FlatVector> stokes{
    {_3,CD::FlatVector{10}}, CD::FlatVector{5} };

  { // check the make-functions
    auto stokes2 = CD::makeDescriptor(
      CD::makeDescriptor(
        CD::makeUniformDescriptor(10, CD::Value{}),
        CD::makeUniformDescriptor(10, CD::Value{}),
        CD::makeUniformDescriptor(10, CD::Value{})
      ),
      CD::makeUniformDescriptor(5, CD::Value{})
    );

    static_assert(std::is_same_v<decltype(stokes),decltype(stokes2)>);
  }

  CD::Array<CD::FlatVector,3> const& velocity = stokes[_0];
  test.check(velocity[0].size() == 10, "v[0].size == 10");
  test.check(velocity[1].size() == 10, "v[1].size == 10");
  test.check(velocity[2].size() == 10, "v[2].size == 10");

  CD::FlatVector const& pressure = stokes[_1];
  test.check(pressure.size() == 5, "p.size == 5");
}

void checkConstructors (Dune::TestSuite& test)
{
  CD::Unknown cd0{};
  CD::Value cd1{};

  // the children might have different types
  CD::Tuple<CD::Unknown,CD::Value> cd2a;
  CD::Tuple<CD::Unknown,CD::Value> cd2b{cd0,cd1};
  CD::Tuple cd2c{cd0,cd1};
  auto cd2d = CD::makeDescriptor(cd0,cd1);
  static_assert(cd2a.size() == 2, "cd2a.size() == 2");
  static_assert(cd2b.size() == 2, "cd2b.size() == 2");
  static_assert(cd2c.size() == 2, "cd2c.size() == 2");
  static_assert(cd2d.size() == 2, "cd2d.size() == 2");
  static_assert(std::is_same_v<decltype(cd2c), decltype(cd2a)>);
  static_assert(std::is_same_v<decltype(cd2d), decltype(cd2a)>);


  // all childresn have the same type, the number of children is static
  CD::Array<CD::Value, 3> cd3a;
  CD::Array<CD::Value, 3> cd3b{cd1};
  CD::Array cd3c{std::integral_constant<std::size_t,3>{},cd1};
  CD::Array cd3d{cd1,cd1,cd1};
  auto cd3e = CD::makeDescriptor(cd1,cd1,cd1);
  static_assert(cd3a.size() == 3, "cd3a.size() == 3");
  static_assert(cd3b.size() == 3, "cd3b.size() == 3");
  static_assert(cd3c.size() == 3, "cd3c.size() == 3");
  static_assert(cd3d.size() == 3, "cd3d.size() == 3");
  static_assert(cd3e.size() == 3, "cd3e.size() == 3");
  static_assert(std::is_same_v<decltype(cd3c), decltype(cd3a)>);
  static_assert(std::is_same_v<decltype(cd3d), decltype(cd3a)>);
  static_assert(std::is_same_v<decltype(cd3e), decltype(cd3a)>);

  // all childresn have the same type, the number of children is a runtime value
  CD::Vector<CD::Value> cd4a(4,cd1);
  CD::Vector<CD::Value> cd4b{cd1,cd1,cd1,cd1};
  CD::Vector cd4c(4,cd1);
  CD::Vector cd4d{cd1,cd1,cd1,cd1};
  test.check(cd4a.size() == 4, "cd4a.size() == 4");
  test.check(cd4b.size() == 4, "cd4b.size() == 4");
  test.check(cd4c.size() == 4, "cd4c.size() == 4");
  test.check(cd4d.size() == 4, "cd4d.size() == 4");
  static_assert(std::is_same_v<decltype(cd4c), decltype(cd4a)>);
  static_assert(std::is_same_v<decltype(cd4d), decltype(cd4a)>);

  // all children are identical, the number of children is static
  // only a single child is stored
  CD::UniformArray<CD::Value, 5> cd5a;
  CD::UniformArray<CD::Value, 5> cd5b(cd1);
  CD::UniformArray cd5c{std::integral_constant<std::size_t,5>{},cd1};
  auto cd5d = makeUniformDescriptor(std::integral_constant<std::size_t,5>{},cd1);
  static_assert(cd5a.size() == 5, "cd5a.size() == 5");
  static_assert(cd5b.size() == 5, "cd5b.size() == 5");
  static_assert(cd5c.size() == 5, "cd5c.size() == 5");
  static_assert(cd5d.size() == 5, "cd5d.size() == 5");
  static_assert(std::is_same_v<decltype(cd5c), decltype(cd5a)>);
  static_assert(std::is_same_v<decltype(cd5d), decltype(cd5a)>);

  // shortcut for uniform arrays storing `Value`
  CD::FlatArray<5> cd5e;
  static_assert(cd5e.size() == 5, "cd5e.size() == 5");
  static_assert(std::is_same_v<decltype(cd5e), decltype(cd5a)>);

  // all children are identical, the number of children is a runtime value
  // only a single child is stored
  CD::UniformVector<CD::Value> cd6a(6);
  CD::UniformVector cd6b(6,cd1);
  auto cd6c = makeUniformDescriptor(6,cd1);
  test.check(cd6a.size() == 6, "cd6a.size() == 6");
  test.check(cd6b.size() == 6, "cd6b.size() == 6");
  test.check(cd6c.size() == 6, "cd6c.size() == 6");
  static_assert(std::is_same_v<decltype(cd6b), decltype(cd6a)>);
  static_assert(std::is_same_v<decltype(cd6c), decltype(cd6a)>);

  // shortcut for uniform vectors storing `Value`
  CD::FlatVector cd6d(6);
  test.check(cd6d.size() == 6, "cd6d.size() == 6");
  static_assert(std::is_same_v<decltype(cd6d), decltype(cd6a)>);
}

int main (int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;
  checkConstructors(test);
  checkHierarchic(test);

  return test.exit();
}
