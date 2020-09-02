// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/functions/common/polymorphicsmallobject.hh>

#include <array>
#include <iostream>
#include <string>
#include <utility>

class Base
{
public:
  virtual ~Base() = default;
  virtual Base* clone() = 0;
  virtual Base* clone(void*) = 0;
  virtual Base* move(void*) = 0;
  virtual int foo() = 0;
  virtual bool checkMemberAlignment() = 0;
};

class Derived
    : public virtual Base
{
  long double dummy_;
  int i_;
public:
  virtual ~Derived() = default;
  Derived(int i = 0) : i_(i) {}
  Base* clone() final { return new Derived{i_}; }
  Base* clone(void* ptr) final { return new (ptr) Derived{this->i_}; }
  Base* move(void* ptr) final { return new (ptr) Derived{std::move(this->i_)}; }
  int foo() final { return i_; }
  bool checkMemberAlignment() final {
    return (reinterpret_cast<std::uintptr_t>(&dummy_) % alignof(long double)) == 0;
  }
};

bool checkTrue(bool value, std::string error)
{
  if (not(value))
    std::cout << "TEST FAILURE:" << error << std::endl;
  return value;
}

template <class Obj>
bool test()
{
  bool success = true;

  std::array<Obj,2> v{Derived{2}, Derived{2}};
  success &= checkTrue(v[0].get().checkMemberAlignment(), "Invalid alignment of member!");
  success &= checkTrue(v[1].get().checkMemberAlignment(), "Invalid alignment of member!");

  Obj obj(Derived{2});
  success &= checkTrue(obj.get().foo() == 2, "Invalid state using constructor with Derived argument!");

  Obj objCopy(obj);
  success &= checkTrue(objCopy.get().foo() == 2, "Invalid state using copy constructor!");

  Obj objMove(std::move(objCopy));
  success &= checkTrue(objMove.get().foo() == 2, "Invalid state using move constructor!");

  Obj objCopyAssign = obj;
  success &= checkTrue(objCopyAssign.get().foo() == 2, "Invalid state using copy assignment!");

  Obj objMoveAssign = std::move(objCopyAssign);
  success &= checkTrue(objMoveAssign.get().foo() == 2, "Invalid state using move assignment!");

  obj = obj;
  success &= checkTrue(obj.get().foo() == 2, "Invalid state using self assignment!");

  return success;
}

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance(argc, argv);
  bool passed = true;

  std::cout << "Testing PolymorphicSmallObject with no buffer" << std::endl;
  passed &= test<Dune::Functions::PolymorphicSmallObject<Base, 0>>();

  std::cout << "Testing PolymorphicSmallObject with buffer" << std::endl;
  constexpr std::size_t OBJSIZE = sizeof(Derived);
  passed &= test<Dune::Functions::PolymorphicSmallObject<Base, OBJSIZE>>();

  if (passed)
    std::cout << "All tests passed" << std::endl;

  return passed ? 0: 1;
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
