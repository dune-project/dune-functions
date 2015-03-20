// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH

#include <type_traits>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/localfunction_imp.hh>
#include <dune/functions/common/polymorphicsmallobject.hh>
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



/**
 * \brief Class storing differentiable functions using type erasure
 *
 */
template<class Range, class Domain, class LocalContext, template<class> class DerivativeTraits, size_t bufferSize>
class LocalFunction< Range(Domain), LocalContext, DerivativeTraits, bufferSize>
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
  using DerivativeInterface = LocalFunction<DerivativeSignature, LocalContext, DerivativeTraits, bufferSize>;

protected:

  using WrapperIf = Imp::LocalFunctionWrapperInterface<Signature, DerivativeInterface, LocalContext>;

  template<class B>
  using WrapperImp = Imp::LocalFunctionWrapperImplementation<Signature, DerivativeInterface, LocalContext, B>;

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
    f_(Imp::TypeErasureDerived<WrapperIf, WrapperImp, typename std::decay<F>::type>(std::forward<F>(f)))
  {}

  LocalFunction() = default;

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
  friend DerivativeInterface derivative(const LocalFunction& t)
  {
    return t.f_.get().derivative();
  }

  void bind(const LocalContext& context)
  {
    return f_.get().bind(context);
  }

  void unbind()
  {
    return f_.get().unbind();
  }

  const LocalContext& localContext() const
  {
    return f_.get().localContext();
  }

private:
  PolymorphicSmallObject<Imp::TypeErasureBase<WrapperIf>, bufferSize > f_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_LOCAL_FUNCTION_HH
