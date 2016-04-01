// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASISTAGS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASISTAGS_HH

#include <type_traits>
#include <dune/common/concept.hh>

namespace Dune {
namespace Functions {

  namespace Concept {

    struct IndexMergingStrategy
    {
      template<typename T>
      auto require(T&& t) -> decltype(
        registerIndexMergingStrategy(t)
      );
    };

    template<typename T>
    static constexpr bool isIndexMergingStrategy()
    {
      return models<Concept::IndexMergingStrategy,T>();
    }

    template<typename T>
    static constexpr bool isIndexMergingStrategy(T&& t)
    {
      return models<Concept::IndexMergingStrategy,std::decay_t<T>>();
    }

  } // namespace Concept


namespace BasisBuilder {


  //! Base class for index merging strategies to simplify detection
  struct IndexMergingStrategy {};

  void registerIndexMergingStrategy(IndexMergingStrategy);

  //! Lexicographic merging of direct children without blocking.
  struct FlatLexicographic
    : public IndexMergingStrategy
  {};

  //! Interleaved merging of direct children without blocking.
  struct FlatInterleaved
    : public IndexMergingStrategy
  {};

  //! Lexicographic merging of direct children with blocking (i.e. creating one block per direct child).
  struct BlockedLexicographic
    : public IndexMergingStrategy
  {};

  //! Interleaved merging of direct children with blocking (i.e. creating blocks at the leaves containing one leaf per child each).
  struct LeafBlockedInterleaved : public IndexMergingStrategy {};


  //! Creates a lexicographic merging of direct children without blocking.
  constexpr FlatLexicographic flatLexicographic()
  {
    return {};
  }

  //! Creates an interleaved merging of direct children without blocking.
  constexpr FlatInterleaved flatInterleaved()
  {
    return {};
  }

  //! Creates a lexicographic merging of direct children with blocking (i.e. creating one block per direct child).
  constexpr BlockedLexicographic blockedLexicographic()
  {
    return {};
  }

  //! Creates an interleaved merging of direct children with blocking (i.e. creating blocks at the leaves containing one leaf per child each).
  constexpr LeafBlockedInterleaved leafBlockedInterleaved()
  {
    return {};
  }

} // end namespace BasisBuilder
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASISTAGS_HH
