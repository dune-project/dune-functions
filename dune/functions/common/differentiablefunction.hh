// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH

#include <type_traits>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/concept.hh>
#include <dune/functions/common/signature.hh>
#include <dune/functions/common/typeerasure.hh>

namespace Dune {
namespace Functions {


namespace Concept {



// Function concept ############################################################
template<class Signature>
struct Function;

template<class Range, class Domain>
struct Function<Range(Domain)>
{
  template<class F>
  auto require(F&& f) -> decltype(
    Range(f(std::declval<Domain>()))
  );
};

/// Check if F models the Function concept with given signature
template<class F, class... T>
static constexpr bool isFunction()
{ return Concept::models<Concept::Function<T...>, F>(); }

/// Check if f models the Function concept with given signature
template<class F, class Signature, template<class> class DerivativeTraits>
static constexpr bool isFunction(F&& f, SignatureTag<Signature, DerivativeTraits>)
{ return Concept::models<Concept::Function<Signature>, F>(); }



// DifferentiableFunction concept ##############################################
template<class Signature, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct DifferentiableFunction;

template<class Range, class Domain, template<class> class DerivativeTraits>
struct DifferentiableFunction<Range(Domain), DerivativeTraits>
{
  using DerivativeSignature = typename SignatureTraits<Range(Domain)>::template DerivativeSignature<DerivativeTraits>;

  template<class F>
  auto require(F&& f) -> decltype(
    Range(f(std::declval<Domain>())),
    requireTrue<isFunction<decltype(derivative(f)), DerivativeSignature>()>()
  );
};

/// Check if F models the DifferentiableFunction concept with given signature
template<class F, class... T>
static constexpr bool isDifferentiableFunction()
{ return Concept::models<Concept::DifferentiableFunction<T...>, F>(); }

/// Check if f models the DifferentiableFunction concept with given signature
template<class F, class Signature, template<class> class DerivativeTraits>
static constexpr bool isDifferentiableFunction(F&& f, SignatureTag<Signature, DerivativeTraits>)
{ return Concept::models<Concept::DifferentiableFunction<Signature, DerivativeTraits>, F>(); }



} // namespace Dune::Functions::Concept




/**
 * Default implementation is empty
 * The actual implementation is only given if Signature is an type
 * describing a function signature as Range(Domain).
 */
template<class Signature, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=56>
class DifferentiableFunction
{};



namespace Imp
{

  /// Traits class providing type information for DifferentiableFunction
  template<class S, template<class> class DerivativeTraits, size_t bufferSize>
  struct DifferentiableFunctionTraits
  {
    /// Signature type
    using Signature = S;

    /// Raw signature with unqualified types
    using RawSignature = typename SignatureTraits<Signature>::RawSignature;

    /// Range type
    using Range = typename SignatureTraits<Signature>::Range;

    /// Domain type
    using Domain = typename SignatureTraits<Signature>::Domain;

    /// Signature of the derivative
    using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

    /// Interface type of the derivative
    using DerivativeInterface = DifferentiableFunction<DerivativeSignature, DerivativeTraits, bufferSize>;

    /// Internal concept type for type erasure
    using Concept = DifferentiableFunctionWrapperInterface<Signature, DerivativeInterface>;

    /// Internal model template for type erasure
    template<class B>
    using Model = DifferentiableFunctionWrapperImplementation<Signature, DerivativeInterface, B>;
  };
}



/**
 * \brief Class storing differentiable functions using type erasure
 *
 */
template<class Range, class Domain, template<class> class DerivativeTraits, size_t bufferSize>
class DifferentiableFunction< Range(Domain), DerivativeTraits, bufferSize> :
  public TypeErasureBase<
    typename Imp::DifferentiableFunctionTraits<Range(Domain), DerivativeTraits, bufferSize>::Concept,
    Imp::DifferentiableFunctionTraits<Range(Domain), DerivativeTraits, bufferSize>::template Model>
{
  using Traits = Imp::DifferentiableFunctionTraits<Range(Domain), DerivativeTraits, bufferSize>;

  using Base = TypeErasureBase<typename Traits::Concept, Traits::template Model>;

  using DerivativeInterface = typename Traits::DerivativeInterface;

public:

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
    Base(std::forward<F>(f))
  {}

  DifferentiableFunction() = default;

  /**
   * \brief Evaluation of wrapped function
   */
  Range operator() (const Domain& x) const
  {
    return this->asInterface().operator()(x);
  }

  /**
   * \brief Get derivative of wrapped function
   *
   * This is a free function that will be found by ADL.
   */
  friend DerivativeInterface derivative(const DifferentiableFunction& t)
  {
    return t.asInterface().derivative();
  }
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
