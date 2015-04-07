// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_SIGNATURE_HH
#define DUNE_FUNCTIONS_COMMON_SIGNATURE_HH

#include <type_traits>

#include <dune/functions/common/defaultderivativetraits.hh>

namespace Dune {
namespace Functions {

/**
 * \brief Helper class to check that F is callable
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
 */
template<class Signature>
struct SignatureTraits;

#ifndef DOXYGEN
/** \brief deduce the signature of the operator() of a class T */
template<class T>
struct SignatureTraits
    : public SignatureTraits<decltype(&T::operator())>
{};

/** \brief deduce the signature of an arbitrary const member function of class C */
template <typename C, typename R, typename D>
struct SignatureTraits<R(C::*)(D) const>
    : public SignatureTraits<R(D)>
{};

/** \brief deduce the signature of an arbitrary member function of class C */
template <typename C, typename R, typename D>
struct SignatureTraits<R(C::*)(D)>
    : public SignatureTraits<R(D)>
{};

/** \brief extract domain and range from a signature (works only for free functions) */
template<class R, class D>
struct SignatureTraits<R(D)>
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


} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_SIGNATURE_HH
