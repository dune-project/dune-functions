// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH

#include <type_traits>

#include <dune/common/typeutilities.hh>

#include <dune/functions/common/typeerasure.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/differentiablefunction.hh>
#include <dune/functions/common/localfunction.hh>
#include <dune/functions/common/functionconcepts.hh>
#include <dune/functions/gridfunctions/localderivativetraits.hh>
#include <dune/functions/gridfunctions/gridfunction_imp.hh>



namespace Dune {
namespace Functions {



/*
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

    /// LocalFunctionTraits associated with this type
    using LocalFunctionTraits = typename Dune::Functions::Imp::LocalFunctionTraits<LocalSignature, Element, LocalDerivativeTraits, bufferSize>;

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
 * \brief Wrapper class for functions defined on a Grid
 *
 * \ingroup FunctionInterface
 *
 * Being defined on a grid means in particular that you can evaluate the function
 * in local coordinates of a given entities of the grid. The set of the entities
 * this function is defined on is given by an EntitySet.
 *
 * This models the \ref Concept::GridFunction<Range(Domain), EntitySet, DerivativeTraits> concept.
 */
template<class Range, class Domain, class ES, template<class> class DerivativeTraits, size_t bufferSize>
class GridFunction<Range(Domain), ES, DerivativeTraits, bufferSize> :
  public TypeErasureBase<
    typename Imp::GridFunctionTraits<Range(Domain), ES, DerivativeTraits, bufferSize>::Concept,
    Imp::GridFunctionTraits<Range(Domain), ES, DerivativeTraits, bufferSize>::template Model>
{
  using Traits = Imp::GridFunctionTraits<Range(Domain), ES, DerivativeTraits, bufferSize>;

  using Base = TypeErasureBase<typename Traits::Concept, Traits::template Model>;

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
    Base(std::forward<F>(f))
  {
    static_assert(Dune::Functions::Concept::isGridFunction<F, Range(Domain), EntitySet, DerivativeTraits>(), "Trying to construct a GridFunction from type that does not model the GridFunction concept");
  }

  GridFunction() = default;

  /**
   * \brief Evaluation of wrapped function
   */
  Range operator() (const Domain& x) const
  {
    return this->asInterface().operator()(x);
  }

  /**
   * \copydoc DifferentiableFunction::derivative
   */
  friend DerivativeInterface derivative(const GridFunction& t)
  {
    return t.asInterface().derivative();
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
    return t.asInterface().wrappedLocalFunction();
  }

  /**
   * \brief Get associated EntitySet
   *
   * This is free function will be found by ADL.
   */
  const EntitySet& entitySet() const
  {
    return this->asInterface().wrappedEntitySet();
  }
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRID_FUNCTION_HH
