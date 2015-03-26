// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_SIGNATURE_HH
#define DUNE_FUNCTIONS_COMMON_SIGNATURE_HH

#include <type_traits>

#include <dune/functions/common/defaultderivativetraits.hh>

namespace Dune {
namespace Functions {



template<class Signature>
struct SignatureTraits;

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
