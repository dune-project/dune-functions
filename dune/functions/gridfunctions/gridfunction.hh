// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH

#include <type_traits>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/localfunction.hh>
#include <dune/functions/common/polymorphicsmallobject.hh>
#include <dune/functions/gridfunctions/gridfunction_imp.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>



namespace Dune {
namespace Functions {



/**
 * Default implementation is empty
 * The actual implementation is only given if Signature is an type
 * describing a function signature as Range(Domain).
 */
template<class Signature, class EntitySet, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=64>
class GridFunction
{};



/**
 * \brief Class storing differentiable functions using type erasure
 *
 */
template<class Range, class Domain, class ES, template<class> class DerivativeTraits, size_t bufferSize>
class GridFunction<Range(Domain), ES, DerivativeTraits, bufferSize>
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
  using DerivativeInterface = GridFunction<DerivativeSignature, ES, DerivativeTraits, bufferSize>;

  using EntitySet = ES;

  // \todo Should this be called local context?
  using Element = typename EntitySet::Element;

  using LocalDomain = typename EntitySet::LocalCoordinate;

  /**
   * \brief Wrapper type of returned local functions
   */
  using LocalFunctionInterface = LocalFunction<Signature, Element, DerivativeTraits, bufferSize>;


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
    f_(Imp::GridFunctionWrapper<Signature, DerivativeInterface, LocalFunctionInterface, EntitySet, typename std::decay<F>::type>(std::forward<F>(f)))
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
  PolymorphicSmallObject<Imp::GridFunctionWrapperBase<Signature, DerivativeInterface, LocalFunctionInterface, EntitySet>, bufferSize > f_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH
