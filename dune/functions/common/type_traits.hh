// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH
#define DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/common/hybridutilities.hh>
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
 * This is useful to avoid failing forwarding constructors
 * for base classes or members.
 */
template<class T, class... Args>
using enableIfConstructible = std::enable_if_t<
  std::is_constructible<T, Args...>::value, int>;



/**
 * \brief Check if type is a statically sized container
 *
 * \ingroup Utility
 *
 * Derives from std::true_type or std::false_type
 */
template<class T>
struct HasStaticSize :
  public IsIntegralConstant<decltype(Dune::Hybrid::size(std::declval<T>()))>
{};

//! A variable template representing the value of \ref HasStaticSize
template<class T>
inline constexpr bool HasStaticSize_v = HasStaticSize<T>::value;


/**
 * \brief Obtain size of statically sized container, or 0 if dynamic size
 *
 * \ingroup Utility
 *
 * Derives from std::integral_constant<std::size_t, size>
 */
template<class T>
struct StaticSizeOrZero :
  public std::conditional_t<HasStaticSize_v<T>,
    decltype(Dune::Hybrid::size(std::declval<T>())),
    std::integral_constant<std::size_t,0>>
{};

/**
 * \brief Obtain size of statically sized container as integral_constant, or fail.
 * \ingroup Utility
 */
template<class T>
using StaticSize = std::enable_if_t<HasStaticSize_v<T>,
  decltype(Dune::Hybrid::size(std::declval<T>()))>;


}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_TYPE_TRAITS_HH
