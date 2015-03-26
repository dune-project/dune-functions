// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
#define DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH


#include <dune/functions/common/concept.hh>
#include <dune/functions/common/signature.hh>

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



}}} // namespace Dune::Functions::Concept

#endif // DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
