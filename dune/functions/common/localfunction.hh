// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH

#include <type_traits>

#include <dune/common/typeutilities.hh>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/localfunction_imp.hh>
#include <dune/functions/common/typeerasure.hh>
#include <dune/functions/common/functionconcepts.hh>



namespace Dune {
namespace Functions {



/*
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
 * \brief Class storing local functions using type erasure
 *
 * \ingroup FunctionInterface
 *
 * \tparam Range Range type
 * \tparam Domain Domain type
 * \tparam LocalContext Type of local context where this function is defined on
 * \tparam DerivativeTraits Traits class to determine range of derivative.
 * \tparam bufferSize Size of stack buffer for small object optimization (defaults to 56)
 *
 *
 * This models the \ref Concept::LocalFunction<Range(Domain), LocalContext, DerivativeTraits> concept.
 * Objects of this type are returned as local functions
 * by the GridFunction wrapper. Notice that the DerivativeTraits type
 * used here should normally be LocalDerivativeTraits<E,GDE> where GDE
 * is the DerivativeTraits type of the corresponding global function.
 * Small object optimization is used to store the given function.
 * If its size exceed \p bufferSize, memory will be allocated dynamically.
 */
template<class Range, class Domain, class LocalContext, template<class> class DerivativeTraits, size_t bufferSize>
class LocalFunction< Range(Domain), LocalContext, DerivativeTraits, bufferSize> :
  public TypeErasureBase<
    typename Imp::LocalFunctionTraits<Range(Domain), LocalContext, DerivativeTraits, bufferSize>::Concept,
    Imp::LocalFunctionTraits<Range(Domain), LocalContext, DerivativeTraits, bufferSize>::template Model>
{
  using Traits = Imp::LocalFunctionTraits<Range(Domain), LocalContext, DerivativeTraits, bufferSize>;

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
  template<class F, disableCopyMove<LocalFunction, F> = 0 >
  LocalFunction(F&& f) :
    Base(std::forward<F>(f))
  {
    static_assert(Dune::Functions::Concept::isLocalFunction<F, Range(Domain), LocalContext>(), "Trying to construct a LocalFunction from type that does not model the LocalFunction concept");
  }

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
   * \ingroup FunctionInterface
   *
   * This is free function will be found by ADL.
   */
  friend DerivativeInterface derivative(const LocalFunction& t)
  {
    return t.asInterface().derivative();
  }

  /**
   * \brief Bind function to a local context
   *
   * You must bind a LocalFunction to a local
   * context before you can evaluate it.
   */
  void bind(const LocalContext& context)
  {
    this->asInterface().bind(context);
  }

  /**
   * \brief Unbind from local context
   */
  void unbind()
  {
    this->asInterface().unbind();
  }

  /** \brief Return if the local function is bound to a grid element
   */
  bool bound() const
  {
    return this->asInterface().bound();
  }

  /**
   * \brief Obtain local contex this LocalFunction is bound to
   */
  const LocalContext& localContext() const
  {
    return this->asInterface().localContext();
  }
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH
