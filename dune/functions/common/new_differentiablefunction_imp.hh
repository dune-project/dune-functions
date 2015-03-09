// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH

#include <string>

#include "concept.hh"

namespace Dune {
namespace Functions {
namespace Imp {


template<class Signature, class DerivativeInterface>
class DifferentiableFunctionWrapperBase
{};

template<class Range, class Domain, class DI>
class DifferentiableFunctionWrapperBase<Range(Domain), DI>
{
public:
  using DerivativeInterface=DI;

  virtual ~DifferentiableFunctionWrapperBase()
  {}

  virtual Range operator() (const Domain& x) const = 0;

  virtual DerivativeInterface wrappedDerivative() const = 0;

  virtual DifferentiableFunctionWrapperBase* clone() const = 0;

  virtual DifferentiableFunctionWrapperBase* clone(void* buffer) const = 0;

  virtual DifferentiableFunctionWrapperBase* move(void* buffer) = 0;
};



template<class Signature, class DerivativeInterface, class FImp>
class DifferentiableFunctionWrapper
{};

template<class Range, class Domain, class DerivativeInterface, class FImp>
class DifferentiableFunctionWrapper< Range(Domain), DerivativeInterface, FImp> :
    public DifferentiableFunctionWrapperBase<Range(Domain), DerivativeInterface>
{
public:

  template<class F>
  DifferentiableFunctionWrapper(F&& f) :
    f_(std::forward<F>(f))
  {}

  virtual Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  virtual DerivativeInterface wrappedDerivative() const
  {
    return DerivativeInterface(derivative(f_));
  }

  virtual DifferentiableFunctionWrapper* clone() const
  {
    return new DifferentiableFunctionWrapper(f_);
  }

  virtual DifferentiableFunctionWrapper* clone(void* buffer) const
  {
    return new (buffer) DifferentiableFunctionWrapper(f_);
  }

  virtual DifferentiableFunctionWrapper* move(void* buffer)
  {
    return new (buffer) DifferentiableFunctionWrapper(std::move(f_));
  }

private:
  FImp f_;
};


template<class Signature, class DerivativeInterface, class FImp>
class NonDifferentiableFunctionWrapper
{};

template<class Range, class Domain, class DerivativeInterface, class FImp>
class NonDifferentiableFunctionWrapper< Range(Domain), DerivativeInterface, FImp> :
  public DifferentiableFunctionWrapperBase<Range(Domain), DerivativeInterface>
{
public:

  template<class F>
  NonDifferentiableFunctionWrapper(F&& f) :
    f_(std::forward<F>(f))
  {}

  virtual Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  virtual DerivativeInterface wrappedDerivative() const
  {
    throw std::string("Derivative not implemented");
  }

  virtual NonDifferentiableFunctionWrapper* clone() const
  {
    return new NonDifferentiableFunctionWrapper(f_);
  }

  virtual NonDifferentiableFunctionWrapper* clone(void* buffer) const
  {
    return new (buffer) NonDifferentiableFunctionWrapper(f_);
  }

  virtual NonDifferentiableFunctionWrapper* move(void* buffer)
  {
    return new (buffer) NonDifferentiableFunctionWrapper(std::move(f_));
  }

private:
  FImp f_;
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
