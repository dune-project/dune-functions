// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/functions/common/concept.hh>
#include <dune/functions/common/test/testsuite.hh>


using namespace Dune::Functions::Concept;



struct HasFoo
{
  template<class T>
  auto require(const T& t) -> decltype(
    t.foo()
  );
};

struct HasBar
{
  template<class T>
  auto require(const T& t) -> decltype(
    t.bar()
  );
};

struct HasFooAndBar1 : Refines<HasFoo>
{
  template<class T>
  auto require(const T& t) -> decltype(
    t.bar()
  );
};

struct HasFooAndBar2 : Refines<HasBar>
{
  template<class T>
  auto require(const T& t) -> decltype(
    t.foo()
  );
};

struct HasFooAndBar3
{
  template<class T>
  auto require(const T& t) -> decltype(
    t.foo(),
    t.bar()
  );
};


template<class T>
struct Foo
{
    T foo() const
    { return T(); }
};

template<class T>
struct Bar
{
    T bar() const
    { return T(); }
};

template<class T>
struct FooBar
{
    T foo() const
    { return T(); }

    T bar() const
    { return T(); }
};




int main ( int argc, char **argv )
try
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  Dune::Functions::TestSuite test;

  test.check(models<HasFoo, Foo<int>>())
      << "models<HasFoo, Foo<int>>() gives wrong result";

  test.check(not models<HasFoo, Bar<int>>())
      << "models<HasFoo, Bar<int>>() gives wrong result";

  test.check(models<HasFoo, FooBar<int>>())
      << "models<HasFoo, FooBar<int>>() gives wrong result";



  test.check(not models<HasBar, Foo<int>>())
      << "models<HasBar, Foo<int>>() gives wrong result";

  test.check(models<HasBar, Bar<int>>())
      << "models<HasBar, Bar<int>>() gives wrong result";

  test.check(models<HasBar, FooBar<int>>())
      << "models<HasBar, FooBar<int>>() gives wrong result";



  test.check(not models<HasFooAndBar1, Foo<int>>())
      << "models<HasFooAndBar1, Foo<int>>() gives wrong result";

  test.check(not models<HasFooAndBar1, Bar<int>>())
      << "models<HasFooAndBar1, Bar<int>>() gives wrong result";

  test.check(models<HasFooAndBar1, FooBar<int>>())
      << "models<HasFooAndBar1, FooBar<int>>() gives wrong result";



  test.check(not models<HasFooAndBar2, Foo<int>>())
      << "models<HasFooAndBar2, Foo<int>>() gives wrong result";

  test.check(not models<HasFooAndBar2, Bar<int>>())
      << "models<HasFooAndBar2, Bar<int>>() gives wrong result";

  test.check(models<HasFooAndBar2, FooBar<int>>())
      << "models<HasFooAndBar2, FooBar<int>>() gives wrong result";



  test.check(not models<HasFooAndBar3, Foo<int>>())
      << "models<HasFooAndBar3, Foo<int>>() gives wrong result";

  test.check(not models<HasFooAndBar3, Bar<int>>())
      << "models<HasFooAndBar3, Bar<int>>() gives wrong result";

  test.check(models<HasFooAndBar3, FooBar<int>>())
      << "models<HasFooAndBar3, FooBar<int>>() gives wrong result";



  return test.exit();
}
catch( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
