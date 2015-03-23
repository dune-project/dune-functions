// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH
#define DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH

#include <type_traits>


namespace Dune {
namespace Functions {


// Helper typedef to remove constructor with forwarding reference from
// overload set for copy and move constructor or assignment.
template<class This, class T>
using disableCopyMove = typename std::enable_if<
  (not(std::is_same<This,  typename std::decay<T>::type>::value)
  and not(std::is_base_of<This,T>::value)), int>::type;



}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH
