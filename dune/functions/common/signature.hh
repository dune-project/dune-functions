// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_SIGNATURE_HH
#define DUNE_FUNCTIONS_COMMON_SIGNATURE_HH

#include <type_traits>
#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/apply.hh>
#include <dune/common/std/type_traits.hh>

#include <dune/functions/common/defaultderivativetraits.hh>

namespace Dune {
namespace Functions {

/**
 * \brief Helper class to check that F is callable
 *
 * \ingroup FunctionUtility
 */
template<typename F>
struct IsCallable;

#ifndef DOXYGEN
template<typename F>
struct IsCallable
{
    struct yes { std::size_t dummy[2]; };
    struct no  { std::size_t dummy[1]; };

    template<typename C>
    static yes test(const decltype(&C::operator()) *);
    template<typename C>
    static no  test(...);

    enum { value = (sizeof(test<F>(0)) == sizeof(yes)) };
};

template<typename R, typename D>
struct IsCallable<R(D)>
{
    enum { value = true };
};

template<typename R, typename D>
struct IsCallable<R(*)(D)>
{
    enum { value = true };
};
#endif

/**
 * \brief Helper class to deduce the signature of a callable
 *
 * \ingroup FunctionUtility
 */
template<class Signature, bool isCallable = IsCallable<Signature>::value >
struct SignatureTraits {};

#ifndef DOXYGEN
/** \brief deduce the signature of the operator() of a class T */
template<class T>
struct SignatureTraits<T, true>
    : public SignatureTraits<decltype(&T::operator()), true>
{};

/** \brief deduce the signature of an arbitrary const member function of class C */
template <typename C, typename R, typename D>
struct SignatureTraits<R(C::*)(D) const, true>
    : public SignatureTraits<R(D), true>
{};

/** \brief deduce the signature of an arbitrary member function of class C */
template <typename C, typename R, typename D>
struct SignatureTraits<R(C::*)(D), true>
    : public SignatureTraits<R(D), true>
{};

/** \brief extract domain and range from a free functions pointer */
template <typename R, typename D>
struct SignatureTraits<R(*)(D), true>
    : public SignatureTraits<R(D), true>
{};

/** \brief extract domain and range from a signature (works only for free functions) */
template<class R, class D>
struct SignatureTraits<R(D), true>
{
    using Range = R;
    using Domain = D;

    using RawRange = typename std::decay<Range>::type;
    using RawDomain = typename std::decay<Domain>::type;

    using RawSignature = RawRange(RawDomain);

    template<template<class> class DerivativeTraits=DefaultDerivativeTraits>
    using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);
};
#endif


template<class Signature, template<class> class DerivativeTraits=DefaultDerivativeTraits>
struct SignatureTag;

/**
 * \brief Tag-class to encapsulate signature information
 *
 * \ingroup FunctionUtility
 *
 * \tparam Range range type
 * \tparam Domain domain type
 * \tparam DerivativeTraits traits template used to determine derivative traits
 */
template<class Range, class Domain, template<class> class DerivativeTraitsT>
struct SignatureTag<Range(Domain), DerivativeTraitsT>
{
  using Signature = Range(Domain);

  template<class T>
  using DerivativeTraits = DerivativeTraitsT<T>;
};



/**
 * \brief Construct SignatureTag for derivative
 *
 * \ingroup FunctionUtility
 *
 * \param tag SignatureTag for a function
 * \returns SignatureTags of the derivative
 */
template<class Range, class Domain, template<class> class DerivativeTraits>
auto derivativeSignatureTag(SignatureTag<Range(Domain), DerivativeTraits> tag)
{
  using DerivativeRange = typename DerivativeTraits<Range(Domain)>::Range;
  return SignatureTag<DerivativeRange(Domain), DerivativeTraits>();
}



/**
 * \brief Construct SignatureTags for derivatives
 *
 * \ingroup FunctionUtility
 *
 * \tparam maxOrder Maximal order of derivatives
 * \param tag SignatureTag for a function
 *
 * \returns Tuple of SignatureTags
 *
 * This constructs an std::tuple of SignatureTags for
 * all derivatives of order 0 up to maxOrder.
 */
template<std::size_t maxOrder, class Signature, template<class> class DerivativeTraits>
auto derivativeSignatureTags(Dune::Functions::SignatureTag<Signature, DerivativeTraits> tag)
{
//  using namespace Dune::Hybrid;
//  using namespace Dune::Std;

  return Hybrid::ifElse(Std::bool_constant<maxOrder==0>(), [&](auto id) {
    // If maxOrder== 0 we just need the given SignatureTag
    return std::make_tuple(tag);
  }, [&](auto id) {
    // else we first construct the tail tuple with SignatureTags for derivatives
    // of order 1 to maxOrder
    auto tailTagsTuple = derivativeSignatureTags<decltype(id,std::size_t(0))(maxOrder-1)>(derivativeSignatureTag(tag));
    // and prepend this with the given SignatureTag.
    // This is done by unpacking the tail tuple with apply().
    return Std::apply([&](auto&&... tailTags){
      return std::make_tuple(tag, tailTags...);
    }, tailTagsTuple);
  });
}



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_SIGNATURE_HH
