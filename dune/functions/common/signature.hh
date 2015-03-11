// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_SIGNATURE_HH
#define DUNE_FUNCTIONS_COMMON_SIGNATURE_HH

#include <type_traits>

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
};



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_SIGNATURE_HH
