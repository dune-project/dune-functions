// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH

#include <dune/common/exceptions.hh>
#include <dune/common/concept.hh>

#include <dune/functions/common/type_traits.hh>


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
    models< HasFreeDerivative, F>() , int>::type = 0>
auto derivativeIfImplemented(const F& f) -> decltype(derivative(f))
{
  return derivative(f);
}



template<class Dummy, class F,
  typename std::enable_if<
    not(models< HasFreeDerivative, F>()) , int>::type = 0>
Dummy derivativeIfImplemented(const F& f)
{
  DUNE_THROW(Dune::NotImplemented, "Derivative not implemented");
}



template<class Signature, class DerivativeInterface>
class DifferentiableFunctionWrapperInterface
{};

// Interface of type erasure wrapper
//
// Notice that the basic interface of polymorphic classes (destructor, clone, ...)
// will be added by the type erasure foundation classes.
template<class Range, class Domain, class DerivativeInterface>
class DifferentiableFunctionWrapperInterface<Range(Domain), DerivativeInterface>
{
public:
  virtual Range operator() (const Domain& x) const = 0;

  virtual DerivativeInterface derivative() const = 0;
};



template<class Signature, class DerivativeInterface, class B>
class DifferentiableFunctionWrapperImplementation
{};

// Implementation of type erasure wrapper
template<class Range, class Domain, class DerivativeInterface, class B>
class DifferentiableFunctionWrapperImplementation< Range(Domain), DerivativeInterface, B> :
  public B
{
public:

  using B::B;
  using Wrapped = typename B::Wrapped;

  virtual Range operator() (const Domain& x) const
  {
    return this->get()(x);
  }

  virtual DerivativeInterface derivative() const
  {
    return derivativeIfImplemented<DerivativeInterface, Wrapped>(this->get());
  }
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
