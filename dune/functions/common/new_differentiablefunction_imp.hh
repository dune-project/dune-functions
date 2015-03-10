// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH

#include <dune/common/exceptions.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/interfaces.hh>

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

template<class Range, class Domain, class DerivativeInterface>
class DifferentiableFunctionWrapperBase<Range(Domain), DerivativeInterface> :
  public PolymorphicType<DifferentiableFunctionWrapperBase<Range(Domain), DerivativeInterface> >
{
public:

  virtual Range operator() (const Domain& x) const = 0;

  virtual DerivativeInterface derivative() const = 0;
};



template<class Signature, class DerivativeInterface, class FImp>
class DifferentiableFunctionWrapper
{};

template<class Range, class Domain, class DerivativeInterface, class FImp>
class DifferentiableFunctionWrapper< Range(Domain), DerivativeInterface, FImp> :
  public DifferentiableFunctionWrapperBase<Range(Domain), DerivativeInterface>
{
public:

  template<class F, disableCopyMove<DifferentiableFunctionWrapper, F> = 0>
  DifferentiableFunctionWrapper(F&& f) :
    f_(std::forward<F>(f))
  {}

  virtual Range operator() (const Domain& x) const
  {
    return f_(x);
  }

  virtual DerivativeInterface derivative() const
  {
    return derivativeIfImplemented<DerivativeInterface, FImp>(f_);
  }

  virtual DifferentiableFunctionWrapper* clone() const
  {
    return new DifferentiableFunctionWrapper(*this);
  }

  virtual DifferentiableFunctionWrapper* clone(void* buffer) const
  {
    return new (buffer) DifferentiableFunctionWrapper(*this);
  }

  virtual DifferentiableFunctionWrapper* move(void* buffer)
  {
    return new (buffer) DifferentiableFunctionWrapper(std::move(*this));
  }

private:
  FImp f_;
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
