// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
#define DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH


#include <dune/functions/common/concept.hh>
#include <dune/functions/common/signature.hh>
#include <dune/functions/gridfunctions/localderivativetraits.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune {
namespace Functions {
namespace Concept {



// Callable concept ############################################################
template<class... Args>
struct Callable
{
  template<class F>
  auto require(F&& f) -> decltype(
    f(std::declval<Args>()...)
  );
};

/// Check if F models the Function concept with given signature
template<class F, class... Args>
static constexpr bool isCallable()
{ return Concept::models<Concept::Callable<Args...>, F>(); }

/// Check if f models the Function concept with given signature
template<class F, class... Args>
static constexpr bool isCallable(F&& f, TypeList<Args...>)
{ return Concept::models<Concept::Callable<Args...>, F>(); }



// Function concept ############################################################
template<class Signature>
struct Function;

template<class Range, class Domain>
struct Function<Range(Domain)> : Refines<Callable<Domain> >
{
  template<class F>
  auto require(F&& f) -> decltype(
    // F models Function<Range(Domain)> if the result of F(Domain) is implicitly convertible to Range
//    requireTrue< std::is_convertible<typename std::result_of<F(Domain)>::type, Range>::value >()
    requireConvertible<typename std::result_of<F(Domain)>::type, Range>()
  );
};

/// Check if F models the Function concept with given signature
template<class F, class Signature>
static constexpr bool isFunction()
{ return Concept::models<Concept::Function<Signature>, F>(); }

/// Check if f models the Function concept with given signature
template<class F, class Signature, template<class> class DerivativeTraits>
static constexpr bool isFunction(F&& f, SignatureTag<Signature, DerivativeTraits>)
{ return Concept::models<Concept::Function<Signature>, F>(); }



// DifferentiableFunction concept ##############################################
template<class Signature, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct DifferentiableFunction;

template<class Range, class Domain, template<class> class DerivativeTraits>
struct DifferentiableFunction<Range(Domain), DerivativeTraits> : Refines<Dune::Functions::Concept::Function<Range(Domain)> >
{
  using DerivativeSignature = typename SignatureTraits<Range(Domain)>::template DerivativeSignature<DerivativeTraits>;

  template<class F>
  auto require(F&& f) -> decltype(
    derivative(f),
    requireTrue<isFunction<decltype(derivative(f)), DerivativeSignature>()>()
  );
};

/// Check if F models the DifferentiableFunction concept with given signature
template<class F, class Signature, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isDifferentiableFunction()
{ return Concept::models<Concept::DifferentiableFunction<Signature, DerivativeTraits>, F>(); }

/// Check if f models the DifferentiableFunction concept with given signature
template<class F, class Signature, template<class> class DerivativeTraits>
static constexpr bool isDifferentiableFunction(F&& f, SignatureTag<Signature, DerivativeTraits>)
{ return Concept::models<Concept::DifferentiableFunction<Signature, DerivativeTraits>, F>(); }



// LocalFunction concept ##############################################
template<class Signature, class LocalContext, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct LocalFunction;

template<class Range, class Domain, class LocalContext, template<class> class DerivativeTraits>
struct LocalFunction<Range(Domain), LocalContext, DerivativeTraits> :
    Refines<Dune::Functions::Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> >
{
  template<class F>
  auto require(F&& f) -> decltype(
    f.bind(std::declval<LocalContext>()),
    f.unbind(),
    f.localContext(),
    requireConvertible<decltype(f.localContext()), LocalContext>()
  );
};

/// Check if F models the LocalFunction concept with given signature and local context
template<class F, class Signature, class LocalContext, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isLocalFunction()
{ return Concept::models<Concept::LocalFunction<Signature, LocalContext, DerivativeTraits>, F>(); }



// EntitySet concept ##############################################
struct EntitySet
{
  template<class E>
  auto require(E&& f) -> decltype(
    requireType<typename E::Element>(),
    requireType<typename E::LocalCoordinate>(),
    requireType<typename E::GlobalCoordinate>()
  );
};

/// Check if F models the GridFunction concept with given signature and entity set
template<class E>
static constexpr bool isEntitySet()
{ return Concept::models<Concept::EntitySet, E>(); }



// GridFunction concept ##############################################
template<class Signature, class EntitySet, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct GridFunction;

template<class Range, class Domain, class EntitySet, template<class> class DerivativeTraits>
struct GridFunction<Range(Domain), EntitySet, DerivativeTraits> :
    Refines<Dune::Functions::Concept::DifferentiableFunction<Range(Domain), DerivativeTraits> >
{
  using LocalSignature = Range(typename EntitySet::LocalCoordinate);
  using LocalContext = typename EntitySet::Element;

  template<class R>
  using LocalDerivativeTraits = typename Dune::Functions::LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits<R>;

  template<class F>
  auto require(F&& f) -> decltype(
    localFunction(f),
    f.entitySet(),
    requireTrue<isLocalFunction<decltype(localFunction(f)), LocalSignature, LocalContext, LocalDerivativeTraits>()> (),
    requireTrue<isEntitySet<EntitySet>()>(),
    requireConvertible<decltype(f.entitySet()), EntitySet>(),
    requireConvertible<typename EntitySet::GlobalCoordinate, Domain>()
  );
};

/// Check if F models the GridFunction concept with given signature and entity set
template<class F, class Signature, class EntitySet, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isGridFunction()
{ return Concept::models<Concept::GridFunction<Signature, EntitySet, DerivativeTraits>, F>(); }



// GridViewFunction concept ##############################################
template<class Signature, class GridView, template<class> class DerivativeTraits = DefaultDerivativeTraits>
struct GridViewFunction;

template<class Range, class Domain, class GridView, template<class> class DerivativeTraits>
struct GridViewFunction<Range(Domain), GridView, DerivativeTraits> :
    Refines<Dune::Functions::Concept::GridFunction<Range(Domain), GridViewEntitySet<GridView,0>, DerivativeTraits>>
{
  template<class F>
  auto require(F&& f) -> decltype(
    0 // We don't need to check any further expressions, because a GridViewFunction is just a GridFunction with a special EntitySet
  );
};

/// Check if F models the GridViewFunction concept with given signature
template<class F, class Signature, class GridView, template<class> class DerivativeTraits = DefaultDerivativeTraits>
static constexpr bool isGridViewFunction()
{ return Concept::models<Concept::GridViewFunction<Signature, GridView, DerivativeTraits>, F>(); }



}}} // namespace Dune::Functions::Concept

#endif // DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
