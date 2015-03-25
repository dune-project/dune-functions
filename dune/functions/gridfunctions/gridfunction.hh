// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH

#include <type_traits>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/localfunction.hh>
#include <dune/functions/common/polymorphicsmallobject.hh>
#include <dune/functions/gridfunctions/localderivativetraits.hh>
#include <dune/functions/gridfunctions/gridfunction_imp.hh>



namespace Dune {
namespace Functions {



/**
 * Default implementation is empty
 * The actual implementation is only given if Signature is an type
 * describing a function signature as Range(Domain).
 */
template<class Signature, class EntitySet, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=56>
class GridFunction
{};



namespace Imp
{

  /// Traits class providing type information for DifferentiableFunction
  template<class S, class ES, template<class> class DerivativeTraits, size_t bufferSize>
  struct GridFunctionTraits :
    DifferentiableFunctionTraits<S, DerivativeTraits, bufferSize>
  {
  protected:
    using Base=DifferentiableFunctionTraits<S, DerivativeTraits, bufferSize>;

  public:
    /// EntitySet the GridFunction lives on
    using EntitySet = ES;

    /// Element type of EntitySet
    using Element = typename EntitySet::Element;

    /// Signature of the derivative
    using DerivativeSignature = typename Base::DerivativeSignature;

    /// Interface type of the derivative
    using DerivativeInterface = GridFunction<DerivativeSignature, ES, DerivativeTraits, bufferSize>;

    /// Signature of the derivative
    using LocalSignature = typename Base::Range(typename EntitySet::LocalCoordinate);

    template<class R>
    using LocalDerivativeTraits = typename Dune::Functions::LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits<R>;

    /// Interface type of the local function
    using LocalFunctionInterface = LocalFunction<LocalSignature, Element, LocalDerivativeTraits, bufferSize>;

    /// Internal concept type for type erasure
    using Concept = GridFunctionWrapperInterface<S, DerivativeInterface, LocalFunctionInterface, ES>;

    /// Internal model template for type erasure
    template<class B>
    using Model = GridFunctionWrapperImplementation<S, DerivativeInterface, LocalFunctionInterface, ES, B>;
  };
}



/**
 * \brief Class storing differentiable functions using type erasure
 *
 */
template<class Range, class Domain, class ES, template<class> class DerivativeTraits, size_t bufferSize>
class GridFunction<Range(Domain), ES, DerivativeTraits, bufferSize>
{
  using Traits = Imp::GridFunctionTraits<Range(Domain), ES, DerivativeTraits, bufferSize>;

  using DerivativeInterface = typename Traits::DerivativeInterface;

  using LocalFunctionInterface = typename Traits::LocalFunctionInterface;

  using EntitySet = typename Traits::EntitySet;

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
  template<class F, disableCopyMove<GridFunction, F> = 0 >
  GridFunction(F&& f) :
    f_(Imp::TypeErasureDerived<typename Traits::Concept, Traits::template Model, typename std::decay<F>::type>(std::forward<F>(f)))
  {}

  GridFunction() = default;

  /**
   * \brief Evaluation of wrapped function
   */
  Range operator() (const Domain& x) const
  {
    return f_.get().operator()(x);
  }

  /**
   * \copydoc DifferentiableFunction::derivative
   */
  friend DerivativeInterface derivative(const GridFunction& t)
  {
    return t.f_.get().derivative();
  }

  /**
   * \brief Get local function of wrapped function
   *
   * This is free function will be found by ADL.
   *
   * Notice that the returned LocalFunction can
   * only be used after it has been bound to a
   * proper local context.
   */
  friend LocalFunctionInterface localFunction(const GridFunction& t)
  {
    return t.f_.get().wrappedLocalFunction();
  }

  /**
   * \brief Get associated EntitySet
   *
   * This is free function will be found by ADL.
   */
  const EntitySet& entitySet() const
  {
    return f_.get().wrappedEntitySet();
  }


private:
  PolymorphicSmallObject<Imp::TypeErasureBase<typename Traits::Concept>, bufferSize > f_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH
