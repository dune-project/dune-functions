// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH
#define DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/common/typeutilities.hh>

namespace Dune {
namespace Functions {


/**
 * \brief Helper to constrain forwarding constructors
 *
 * \ingroup Utility
 *
 * Helper typedef to remove constructor with forwarding reference from
 * overload set for type is not constructible from argument list.
 * This is useful to avoid failng forwarding constructors
 * for base classes or members.
 */
template<class T, class... Args>
using enableIfConstructible = typename std::enable_if<
  std::is_constructible<T, Args...>::value, int>::type;



namespace Imp {

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
 * \ingroup Utility
 *
 * Derives from std::true_type or std::false_type
 */
template<class T>
struct HasStaticSize :
  public decltype(Imp::hasStaticSize((typename std::decay<T>::type*)(nullptr), PriorityTag<42>()))
{};



/**
 * \brief Obtain size of statically sized container
 *
 * \ingroup Utility
 *
 * Derives from std::integral_constant<std::size_t, size>
 */
template<class T>
struct StaticSize :
  public decltype(Imp::staticSize((typename std::decay<T>::type*)(nullptr), PriorityTag<42>()))
{};



}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH
