// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH

#include <dune/common/exceptions.hh>

#include "concept.hh"

namespace Dune {
namespace Functions {
namespace Imp {

/**
 * A concept describing types that have a derivative() method found by ADL
 */
struct HasFreeDerivative
{
  template<class F>
  auto require(F&& f) -> decltype(
    derivative(f)
  );
};



template<class Dummy, class F,
  typename std::enable_if<
    Dune::Functions::Concept::models< HasFreeDerivative, F>() , int>::type = 0>
auto derivativeIfImplemented(const F& f) -> decltype(derivative(f))
{
  return derivative(f);
}



template<class Dummy, class F,
  typename std::enable_if<
    not(Dune::Functions::Concept::models< HasFreeDerivative, F>()) , int>::type = 0>
Dummy derivativeIfImplemented(const F& f)
{
  DUNE_THROW(Dune::NotImplemented, "Derivative not implemented");
}



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
    return derivativeIfImplemented<DerivativeInterface, FImp>(f_);
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



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
