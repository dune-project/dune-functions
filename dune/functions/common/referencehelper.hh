// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_REFERENCE_HELPER_HH
#define DUNE_FUNCTIONS_COMMON_REFERENCE_HELPER_HH

#warning The header dune/functions/common/referencehelper.hh is deprecated and will be removed after release 2.9. Include dune/common/referencehelper.hh instead.

#include <type_traits>

#include <dune/common/referencehelper.hh>




namespace Dune {
namespace Functions {


/**
 * \brief This is an alias for Dune::IsReferenceWrapper_v
 * \deprecated Use Dune::IsReferenceWrapper_v instead
 */
template<class T>
[[deprecated("Use Dune::IsReferenceWrapper_v instead. Will be removed after release 2.9.")]]
constexpr bool
IsReferenceWrapper_v = Dune::IsReferenceWrapper_v<T>;


/**
 * \brief This is an alias for Dune::resolveRef
 * \deprecated Use Dune::resolveRef instead
 */
template<class T>
decltype(auto)
resolveRef
[[deprecated("Use Dune::resolveRef instead. Will be removed after release 2.9.")]]
(T&& t)
{
  return Dune::resolveRef(std::forward<T>(t));
}

/**
 * \brief This is an alias for Dune::ResolveRef_t
 * \deprecated Use Dune::resolveRef instead
 */
template<class T>
using ResolveRef_t
[[deprecated("Use Dune::ResolveRef_t instead. Will be removed after release 2.9.")]]
  = Dune::ResolveRef_t<T>;


}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_REFERENCE_HELPER_HH
