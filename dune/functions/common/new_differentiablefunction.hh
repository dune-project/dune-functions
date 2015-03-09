// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH

#include <type_traits>

#include "defaultderivativetraits.hh"
#include "new_differentiablefunction_imp.hh"
#include "smallobject.hh"

#include "concept.hh"

namespace Dune {
namespace Functions {
namespace Concept {

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

/**
 * A concept describing types that have a member function foo.derivative()
 *
 * \todo can we remove this?
 */
struct HasMemberDerivative
{
  template<class F>
  auto require(F&& f) -> decltype(
    f.derivative()
  );
};

} // namespace Concept



/**
 * Default implementation is empty
 * The actual implementation is only given if Signature is an type
 * describing a function signature as Range(Domain).
 */
template<class Signature, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=64>
class DifferentiableFunction
{};



/**
 * \brief Class storing differentiable functions using type erasure
 *
 */
template<class Range, class Domain, template<class> class DerivativeTraits, size_t bufferSize>
class DifferentiableFunction< Range(Domain), DerivativeTraits, bufferSize>
{
  using RawDomain = typename std::decay<Domain>::type;
public:

  /**
   * \brief Signature of wrapped functions
   */
  using Signature = Range(Domain);

  /**
   * \brief Signature of derivative of wrapped functions
   */
  using DerivativeSignature = typename DerivativeTraits<Range(RawDomain)>::Range(Domain);

  /**
   * \brief Wrapper type of returned derivatives
   */
  using DerivativeInterface = DifferentiableFunction<DerivativeSignature, DerivativeTraits, bufferSize>;

  /**
   * \brief Construct from function providing a free derivative() function
   *
   * \tparam F Function type
   *
   * \param f Function of type F
   *
   * This constructor is used if there is a free derivative() method
   * to be found by ADL.
   */
  template<class F,
    typename std::enable_if<
      Dune::Functions::Concept::models< Dune::Functions::Concept::HasFreeDerivative, F>()
      and not(std::is_same<DifferentiableFunction, typename std::decay<F>::type>::value), int>::type = 0>
  DifferentiableFunction(F&& f) :
    f_(Imp::DifferentiableFunctionWrapper<Signature, DerivativeInterface, typename std::decay<F>::type>(std::forward<F>(f)))
  {}

  /**
   * \brief Construct from function not providing a free derivative() function
   *
   * \tparam F Function type
   *
   * \param f Function of type F
   *
   * This constructor is used if there is no free derivative() method
   * to be found by ADL. In this case calling derivative for the wrapper
   * throws an exception.
   */
  template<class F,
    typename std::enable_if< not(Dune::Functions::Concept::models< Dune::Functions::Concept::HasFreeDerivative, F>()), int>::type = 0>
  DifferentiableFunction(F&& f) :
    f_(Imp::NonDifferentiableFunctionWrapper<Signature, DerivativeInterface, typename std::decay<F>::type>(std::forward<F>(f)))
  {}


  DifferentiableFunction() = default;

  DifferentiableFunction(const DifferentiableFunction& other) = default;

  DifferentiableFunction(DifferentiableFunction&& other) = default;

  DifferentiableFunction& operator=(const DifferentiableFunction& other) = default;

  DifferentiableFunction& operator=(DifferentiableFunction&& other) = default;

  /**
   * \brief Evaluation of wrapped function
   */
  Range operator() (const Domain& x) const
  {
    return f_.get().operator()(x);
  }

  /**
   * \brief Get derivative of wrapped function
   *
   * This is free function will be found by ADL.
   */
  friend DerivativeInterface derivative(const DifferentiableFunction& t)
  {
    return t.f_.get().wrappedDerivative();
  }

private:
  SmallObject<Imp::DifferentiableFunctionWrapperBase<Signature, DerivativeInterface>, bufferSize > f_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
