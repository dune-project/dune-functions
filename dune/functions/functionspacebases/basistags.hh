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

  /**
   * \brief Base class for index merging strategies to simplify detection
   *
   * \ingroup FunctionSpaceBasesUtilities
   */
  struct IndexMergingStrategy {};

  void registerIndexMergingStrategy(IndexMergingStrategy);

  /**
   * \brief Lexicographic merging of direct children without blocking.
   *
   * \ingroup FunctionSpaceBasesUtilities
   *
   * Example: For two children with indices
   *
   *   {(0,i0), (0,i1), (1,i2)} and {(0,k0), (1,k1), (2,k2)}
   *
   * the merged indices will be
   *
   *   {(0,i0), (0,i1), (1,i2), (2,k0), (3,k1), (4,k2)}
   */
  struct FlatLexicographic
    : public IndexMergingStrategy
  {};

  /**
   * \brief Interleaved merging of direct children without blocking.
   *
   * \ingroup FunctionSpaceBasesUtilities
   *
   * Example: For two children with indices
   *
   *   {(0,i0), (1,i1), (2,i2)} and {(0,i0), (1,i1), (2,i2)}
   *
   * the merged indices will be
   *
   *   {(0,i0), (2,i1), (4,i2), (1,k0), (3,k1), (5,k2)}
   */
  struct FlatInterleaved
    : public IndexMergingStrategy
  {};

  /**
   * \brief Lexicographic merging of direct children with blocking (i.e. creating one block per direct child).
   *
   * \ingroup FunctionSpaceBasesUtilities
   *
   * Example: For two children with indices
   *
   *   {(0,i0), (0,i1), (1,i2)} and {(0,k0), (1,k1), (2,k2)}
   *
   * the merged indices will be
   *
   *   {(0,0,i0), (0,0,i1), (0,1,i2), (1,0,k0), (1,1,k1), (1,2,k2)}
   */
  struct BlockedLexicographic
    : public IndexMergingStrategy
  {};

  /**
   * \brief Interleaved merging of direct children with blocking (i.e. creating blocks at the leaves containing one leaf per child each).
   *
   * \ingroup FunctionSpaceBasesUtilities
   *
   * Example: For two children with indices
   *
   *   {(0,i0), (1,i1), (2,i2)} and {(0,i0), (1,i1), (2,i2)}
   *
   * the merged indices will be
   *
   *   {(0,i0,0), (1,i1,0), (2,i2,0), (0,i0,1), (1,i1,1), (2,i2,1)}
   */
  struct LeafBlockedInterleaved : public IndexMergingStrategy {};


  /**
   * \brief Creates a lexicographic merging of direct children without blocking.
   *
   * \ingroup FunctionSpaceBasesUtilities
   */
  constexpr FlatLexicographic flatLexicographic()
  {
    return {};
  }

  /**
   * \brief Creates an interleaved merging of direct children without blocking.
   *
   * \ingroup FunctionSpaceBasesUtilities
   */
  constexpr FlatInterleaved flatInterleaved()
  {
    return {};
  }

  /**
   * \brief Creates a lexicographic merging of direct children with blocking (i.e. creating one block per direct child).
   *
   * \ingroup FunctionSpaceBasesUtilities
   */
  constexpr BlockedLexicographic blockedLexicographic()
  {
    return {};
  }

  /**
   * \brief Creates an interleaved merging of direct children with blocking (i.e. creating blocks at the leaves containing one leaf per child each).
   *
   * \ingroup FunctionSpaceBasesUtilities
   */
  constexpr LeafBlockedInterleaved leafBlockedInterleaved()
  {
    return {};
  }

} // end namespace BasisBuilder
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASISTAGS_HH
