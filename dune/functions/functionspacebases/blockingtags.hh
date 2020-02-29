// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BLOCKINGTAGS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BLOCKINGTAGS_HH

#include <tuple>
#include <type_traits>

#include <dune/common/indices.hh>

namespace Dune {

// Lightweight representation of (hierarchic) blocking structure
namespace BlockingTag {

  struct Unknown {};

  //! Blocking structure corresponding to IndexMergingStrategies FlatLexicographic and FlatInterleaved
  struct Flat
  {
    template <class Index>
    constexpr Flat operator[](const Index&) const { return {}; }
  };

  //! Blocking structure corresponding to IndexMergingStrategy BlockedInterleaved
  template <std::size_t N>
  struct LeafBlocked
  {
    template <class Index>
    constexpr Flat operator[](const Index&) const { return {}; }
  };

  //! Blocking structure corresponding to IndexMergingStrategy BlockedLexicographic
  template <class Tag0, class... Tags>
  struct Blocked
  {
    template <std::size_t i>
    constexpr auto operator[](index_constant<i>) const { return std::get<i>(t); }
    constexpr Tag0 operator[](std::size_t /*i*/) const { return {}; }

    std::tuple<Tag0,Tags...> t;
  };


  namespace Impl_
  {
    template <class Tag, class Seq>
    struct PowerBlocked;

    template <class Tag, std::size_t... n>
    struct PowerBlocked<Tag, std::index_sequence<n...>>
    {
      template <std::size_t i>
      using expand = Tag;

      using type = Dune::BlockingTag::Blocked<expand<n>...>;
    };
  }

  // helper type for repeated tags in a blocking
  template <class Tag, std::size_t n>
  using PowerBlocked = typename Impl_::PowerBlocked<Tag, std::make_index_sequence<n>>::type;

}} // end namespace Dune::BlockingTag

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BLOCKINGTAGS_HH
