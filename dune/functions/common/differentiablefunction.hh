// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH

#include <type_traits>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/smallobject.hh>
#include <dune/functions/common/concept.hh>
#include <dune/functions/common/signature.hh>

namespace Dune {
namespace Functions {



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
public:

  /**
   * \brief Signature of wrapped functions
   */
  using Signature = Range(Domain);

  /**
   * \brief Raw signature of wrapped functions without possible const and reference qualifiers
   */
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;

  /**
   * \brief Signature of derivative of wrapped functions
   */
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  /**
   * \brief Wrapper type of returned derivatives
   */
  using DerivativeInterface = DifferentiableFunction<DerivativeSignature, DerivativeTraits, bufferSize>;

  /**
   * \brief Construct from function
   *
   * \tparam F Function type
   *
   * \param f Function of type F
   *
   * Calling derivative(DifferentiableFunction) will result in an exception
   * if the passed function does provide a free derivative() function
   * found via ADL.
   */
  template<class F, disableCopyMove<DifferentiableFunction, F> = 0 >
  DifferentiableFunction(F&& f) :
    f_(Imp::DifferentiableFunctionWrapper<Signature, DerivativeInterface, typename std::decay<F>::type>(std::forward<F>(f)))
  {}

  DifferentiableFunction() = default;

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
    return t.f_.get().derivative();
  }

private:
  SmallObject<Imp::DifferentiableFunctionWrapperBase<Signature, DerivativeInterface>, bufferSize > f_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
