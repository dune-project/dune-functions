// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH

#include <type_traits>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/localfunction_imp.hh>
#include <dune/functions/common/typeerasure.hh>



namespace Dune {
namespace Functions {



/**
 * Default implementation is empty
 * The actual implementation is only given if Signature is an type
 * describing a function signature as Range(Domain).
 */
template<class Signature, class LocalContext, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=56>
class LocalFunction
{};



namespace Imp
{

  /// Traits class providing type information for DifferentiableFunction
  template<class S, class L, template<class> class DerivativeTraits, size_t bufferSize>
  struct LocalFunctionTraits :
    DifferentiableFunctionTraits<S, DerivativeTraits, bufferSize>
  {
  protected:
    using Base=DifferentiableFunctionTraits<S, DerivativeTraits, bufferSize>;

  public:
    /// LocalContext type
    using LocalContext = L;

    /// Signature of the derivative
    using DerivativeSignature = typename Base::DerivativeSignature;

    /// Interface type of the derivative
    using DerivativeInterface = LocalFunction<DerivativeSignature, L, DerivativeTraits, bufferSize>;

    /// Internal concept type for type erasure
    using Concept = LocalFunctionWrapperInterface<S, DerivativeInterface, L>;

    /// Internal model template for type erasure
    template<class B>
    using Model = LocalFunctionWrapperImplementation<S, DerivativeInterface, L, B>;
  };
}



/**
 * \brief Class storing differentiable functions using type erasure
 *
 */
template<class Range, class Domain, class LocalContext, template<class> class DerivativeTraits, size_t bufferSize>
class LocalFunction< Range(Domain), LocalContext, DerivativeTraits, bufferSize> :
  public TypeErasure<
    typename Imp::LocalFunctionTraits<Range(Domain), LocalContext, DerivativeTraits, bufferSize>::Concept,
    Imp::LocalFunctionTraits<Range(Domain), LocalContext, DerivativeTraits, bufferSize>::template Model>
{
  using Traits = Imp::LocalFunctionTraits<Range(Domain), LocalContext, DerivativeTraits, bufferSize>;

  using Base = TypeErasure<typename Traits::Concept, Traits::template Model>;

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
  template<class F, disableCopyMove<LocalFunction, F> = 0 >
  LocalFunction(F&& f) :
    Base(std::forward<F>(f))
  {}

  LocalFunction() = default;

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
   * This is free function will be found by ADL.
   */
  friend DerivativeInterface derivative(const LocalFunction& t)
  {
    return t.asInterface().derivative();
  }

  void bind(const LocalContext& context)
  {
    this->asInterface().bind(context);
  }

  void unbind()
  {
    this->asInterface().unbind();
  }

  const LocalContext& localContext() const
  {
    return this->asInterface().localContext();
  }
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH
