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





namespace Imp {

  // Helper class for multi-stage SFINAE based overload
  // resolution. You can tag overloads by priority.
  // The matching overload with the highest prioriry
  // will be selected. Be carefull to always use
  // a large enough priority value when trying to call
  // the overloaded function.
  template<std::size_t i>
  struct PriorityTag : public PriorityTag<i-1>
  {};

  template<>
  struct PriorityTag<0>
  {};



  // As a last resort try if there's a static constexpr size()
  template<class T>
  constexpr auto staticSize(const T*, const PriorityTag<0>&)
    -> decltype(std::integral_constant<std::size_t,T::size()>())
  {
    return {};
  }

  // Try if class has constexpr default constructor and size method
  template<class T>
  constexpr auto staticSize(const T*, const PriorityTag<1>&)
    -> decltype(std::integral_constant<std::size_t,T().size()>())
  {
    return {};
  }

  // Try if tuple_size is implemented for class
  template<class T>
  constexpr auto staticSize(const T*, const PriorityTag<2>&)
    -> decltype(std::integral_constant<std::size_t,std::tuple_size<T>::value>())
  {
    return {};
  }

  template<class T>
  constexpr std::false_type hasStaticSize(const T* t, const PriorityTag<0>& p)
  {
    return {};
  }

  template<class T>
  constexpr auto hasStaticSize(const T* t, const PriorityTag<1>& p)
    -> decltype(staticSize(t ,PriorityTag<42>()), std::true_type())
  {
    return {};
  }

}



/**
 * \brief Check if type is a statically sized container
 *
 * Derives from std::true_type or std::false_type
 */
template<class T>
struct HasStaticSize :
  public decltype(Imp::hasStaticSize((typename std::decay<T>::type*)(nullptr), Imp::PriorityTag<42>()))
{};



/**
 * \brief Obtain size of statically sized container
 *
 * Derives from std::integral_constant<std::size_t, size>
 */
template<class T>
struct StaticSize :
  public decltype(Imp::staticSize((typename std::decay<T>::type*)(nullptr), Imp::PriorityTag<42>()))
{};



}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH
