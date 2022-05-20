// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
#define DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH

#include <dune/common/typelist.hh>
#include <dune/common/concept.hh>

#include <dune/functions/common/signature.hh>
#include <dune/functions/gridfunctions/localderivativetraits.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune {
namespace Functions {
namespace Concept {

using namespace Dune::Concept;



// Callable concept ############################################################


/**
 * \brief Concept objects that can be called with given argument list
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Args Argument list for function call
 */
template<class... Args>
struct Callable
{
  template<class F>
  auto require(F&& f) -> decltype(
    f(std::declval<Args>()...)
  );
};

/**
 * \brief Check if f is callable with given argument list
 *
 * \ingroup FunctionConcepts
 * \ingroup Utility
 */
template<class F, class... Args>
static constexpr auto isCallable()
{ return models<Concept::Callable<Args...>, F>(); }

/**
 * \brief Check if f is callable with given argument list
 *
 * \ingroup FunctionConcepts
 * \ingroup Utility
 */
template<class F, class... Args>
static constexpr auto isCallable(F&&, Args&&...)
{
  return models<Concept::Callable<Args&&...>, F>();
}



// Function concept ############################################################
template<class Signature>
struct Function;

/**
 * \brief Concept for a function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 */
template<class Range, class Domain>
struct Function<Range(Domain)> : Refines<Callable<Domain> >
{
  template<class F>
  auto require(F&& f) -> decltype(
    // F models Function<Range(Domain)> if the result of F(Domain) is implicitly convertible to Range
    requireConvertible<Range>(f(std::declval<Domain>()))
  );
};

/// Check if F models the Function concept with given signature \ingroup FunctionConcepts
template<class F, class Signature>
static constexpr bool isFunction()
{ return models<Concept::Function<Signature>, F>(); }

/// Check if f models the Function concept with given signature \ingroup FunctionConcepts
template<class F, class Signature, template<class> class DerivativeTraits>
static constexpr bool isFunction(F&& f, SignatureTag<Signature, DerivativeTraits>)
{ return models<Concept::Function<Signature>, F>(); }



// DifferentiableFunction concept ##############################################
template<class Signature, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct DifferentiableFunction;

/**
 * \brief Concept for a differentiable function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * The derivative range is derived from the provided \p DerivativeTraits
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template<class Range, class Domain, template<class> class DerivativeTraits>
struct DifferentiableFunction<Range(Domain), DerivativeTraits> : Refines<Dune::Functions::Concept::Function<Range(Domain)> >
{
  using DerivativeSignature = typename SignatureTraits<Range(Domain)>::template DerivativeSignature<DerivativeTraits>;

  template<class F>
  auto require(F&& f) -> decltype(
    derivative(f),
    requireConcept<Function<DerivativeSignature>>(derivative(f))
  );
};

/// Check if F models the DifferentiableFunction concept with given signature \ingroup FunctionConcepts
template<class F, class Signature, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isDifferentiableFunction()
{ return models<Concept::DifferentiableFunction<Signature, DerivativeTraits>, F>(); }

/// Check if f models the DifferentiableFunction concept with given signature \ingroup FunctionConcepts
template<class F, class Signature, template<class> class DerivativeTraits>
static constexpr bool isDifferentiableFunction(F&& f, SignatureTag<Signature, DerivativeTraits>)
{ return models<Concept::DifferentiableFunction<Signature, DerivativeTraits>, F>(); }



// LocalFunction concept ##############################################
template<class Signature, class LocalContext>
struct LocalFunction;

/**
 * \brief Concept for a local function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam LocalContext The local context this function is defined on
 */
template<class Range, class Domain, class LocalContext>
struct LocalFunction<Range(Domain), LocalContext> :
      Refines<Dune::Functions::Concept::Function<Range(Domain)> >
{
  template<class F>
  auto require(F&& f) -> decltype(
    f.bind(std::declval<LocalContext>()),
    f.unbind(),
    requireConvertible<bool>(f.bound()),
    f.localContext(),
    requireConvertible<LocalContext>(f.localContext())
  );
};

/// Check if F models the LocalFunction concept with given signature and local context \ingroup FunctionConcepts
template<class F, class Signature, class LocalContext>
static constexpr bool isLocalFunction()
{ return models<Concept::LocalFunction<Signature, LocalContext>, F>(); }


// DifferentiableLocalFunction concept ##############################################
template<class Signature, class LocalContext, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct DifferentiableLocalFunction;

/**
 * \brief Concept for a differentiable local function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * The derivative range is derived from the provided \p DerivativeTraits
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam LocalContext The local context this function is defined on
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template<class Range, class Domain, class LocalContext, template<class> class DerivativeTraits>
struct DifferentiableLocalFunction<Range(Domain), LocalContext, DerivativeTraits> :
    Refines<
      Dune::Functions::Concept::DifferentiableFunction<Range(Domain), DerivativeTraits>,
      Dune::Functions::Concept::LocalFunction<Range(Domain),LocalContext>
      >
{
  template<class F>
  auto require(F&& f) -> decltype(
    f.bind(std::declval<LocalContext>()),
    f.unbind(),
    f.localContext(),
    requireConvertible<LocalContext>(f.localContext())
  );
};

/// Check if F models the DifferentiableLocalFunction concept with given signature and local context \ingroup FunctionConcepts
template<class F, class Signature, class LocalContext, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isDifferentiableLocalFunction()
{ return models<Concept::DifferentiableLocalFunction<Signature, LocalContext, DerivativeTraits>, F>(); }


// EntitySet concept ##############################################

/**
 * \brief Concept for an entity set for a \ref Concept::GridFunction<Range(Domain), EntitySet, DerivativeTraits>
 *
 * \ingroup FunctionConcepts
 *
 * This describes the set of entities on which a grid function
 * can be localized.
 *
 */
struct EntitySet
{
  template<class E>
  auto require(E&& f) -> decltype(
    requireType<typename E::Element>(),
    requireType<typename E::LocalCoordinate>(),
    requireType<typename E::GlobalCoordinate>()
  );
};

/// Check if F models the GridFunction concept with given signature and entity set \ingroup FunctionConcepts
template<class E>
static constexpr bool isEntitySet()
{ return models<Concept::EntitySet, E>(); }



// GridFunction concept ##############################################
template<class Signature, class EntitySet>
struct GridFunction;

/**
 * \brief Concept for a grid function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam EntitySet Set of entities on which the function can be localized
 */
template<class Range, class Domain, class EntitySet>
struct GridFunction<Range(Domain), EntitySet> :
  Refines<Dune::Functions::Concept::Function<Range(Domain)> >
{
  using LocalSignature = Range(typename EntitySet::LocalCoordinate);
  using LocalContext = typename EntitySet::Element;

  template<class F>
  auto require(F&& f) -> decltype(
    localFunction(f),
    f.entitySet(),
    requireConcept<LocalFunction<LocalSignature, LocalContext>>(localFunction(f)),
    requireConcept<Concept::EntitySet, EntitySet>(),
    requireConvertible<EntitySet>(f.entitySet()),
    requireConvertible<typename EntitySet::GlobalCoordinate, Domain>()
  );
};

/// Check if F models the GridFunction concept with given signature and entity set \ingroup FunctionConcepts
template<class F, class Signature, class EntitySet>
static constexpr bool isGridFunction()
{ return models<Concept::GridFunction<Signature, EntitySet>, F>(); }


// DifferentiableGridFunction concept ##############################################
template<class Signature, class EntitySet, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct DifferentiableGridFunction;

/**
 * \brief Concept for a differentiable grid function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * The derivative range is derived from the provided \p DerivativeTraits.
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam EntitySet Set of entities on which the function can be localized
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template<class Range, class Domain, class EntitySet, template<class> class DerivativeTraits>
struct DifferentiableGridFunction<Range(Domain), EntitySet, DerivativeTraits> :
  Refines<
    Dune::Functions::Concept::DifferentiableFunction<Range(Domain), DerivativeTraits>,
    Dune::Functions::Concept::GridFunction<Range(Domain),EntitySet>
    >
{
  using LocalSignature = Range(typename EntitySet::LocalCoordinate);
  using LocalContext = typename EntitySet::Element;

  template<class R>
  using LocalDerivativeTraits = typename Dune::Functions::LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits<R>;

  template<class F>
  auto require(F&& f) -> decltype(
    requireConcept<DifferentiableLocalFunction<LocalSignature, LocalContext, LocalDerivativeTraits>>(localFunction(f))
  );
};

/// Check if F models the DifferentiableGridFunction concept with given signature and entity set \ingroup FunctionConcepts
template<class F, class Signature, class EntitySet, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isDifferentiableGridFunction()
{ return models<Concept::DifferentiableGridFunction<Signature, EntitySet, DerivativeTraits>, F>(); }



// GridViewFunction concept ##############################################
template<class Signature, class GridView>
struct GridViewFunction;

/**
 * \brief Concept for a grid view function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * This exactly the \ref Concept::GridFunction<Range(Domain), EntitySet>
 * concept with a \ref GridViewEntitySet as \p EntitySet.
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam GridView GridView on which the function can be localized
 */
template<class Range, class Domain, class GridView>
struct GridViewFunction<Range(Domain), GridView> :
  Refines<Dune::Functions::Concept::GridFunction<Range(Domain), GridViewEntitySet<GridView,0>>>
{
  template<class F>
  auto require(F&& f) -> decltype(
    0 // We don't need to check any further expressions, because a GridViewFunction is just a GridFunction with a special EntitySet
  );
};

/// Check if F models the GridViewFunction concept with given signature \ingroup FunctionConcepts
template<class F, class Signature, class GridView>
static constexpr bool isGridViewFunction()
{ return models<Concept::GridViewFunction<Signature, GridView>, F>(); }


// DifferentiableGridViewFunction concept ##############################################
template<class Signature, class GridView, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct DifferentiableGridViewFunction;

/**
 * \brief Concept for a differentiable grid view function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * This exactly the \ref Concept::GridFunction<Range(Domain), EntitySet, DerivativeTraits>
 * concept with a \ref GridViewEntitySet as \p EntitySet.
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam GridView GridView on which the function can be localized
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template<class Range, class Domain, class GridView, template<class> class DerivativeTraits>
struct DifferentiableGridViewFunction<Range(Domain), GridView, DerivativeTraits> :
  Refines<Dune::Functions::Concept::DifferentiableGridFunction<Range(Domain), GridViewEntitySet<GridView,0>, DerivativeTraits>>
{
  template<class F>
  auto require(F&& f) -> decltype(
    0 // We don't need to check any further expressions, because a GridViewFunction is just a GridFunction with a special EntitySet
  );
};

/// Check if F models the DifferentiableGridViewFunction concept with given signature \ingroup FunctionConcepts
template<class F, class Signature, class GridView, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isDifferentiableGridViewFunction()
{ return models<Concept::DifferentiableGridViewFunction<Signature, GridView, DerivativeTraits>, F>(); }



}}} // namespace Dune::Functions::Concept

#endif // DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
