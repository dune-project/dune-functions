// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TYPELIST_HH
#define DUNE_FUNCTIONS_COMMON_TYPELIST_HH

#include <type_traits>


namespace Dune {
namespace Functions {

  /**
   * \brief A simple type list
   *
   * If T... is not empty it provides access to Head and Tail.
   */
  template<class... T>
  struct TypeList
  {};

  template<class T0, class... T>
  struct TypeList<T0, T...>
  {
    typedef T0 Head;
    typedef TypeList<T...> Tail;
  };


  // Simplify recognition of lists
  template<class T>
  struct IsTypeList : std::false_type {};

  template<class... T>
  struct IsTypeList<TypeList<T...> > : std::true_type {};

  /// Constexpr function which is true for type lists
  template<class T>
  static constexpr bool isTypeList()
  { return IsTypeList<T>(); }

  /// Constexpr function which is true for an empty type lists
  template<class T>
  static constexpr bool isEmptyTypeList()
  { return isTypeList<T>() and std::is_same<T, TypeList<> >(); }

}} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_TYPELIST_HH
