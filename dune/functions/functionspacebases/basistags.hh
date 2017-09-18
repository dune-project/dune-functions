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
   * Example: For two children {f} and {g} with multi-indices (all i*,k* can be multi-indices, themselves)
   *
   * function in {f}| index
   * ---------------|------
   *  f_0           |(0,i0)
   *  f_1           |(0,i1)
   *  f_2           |(1,i2)
   *
   * function in {g}| index
   * ---------------|------
   *  g_0           | (0,k0)
   *  g_1           | (1,k1)
   *  g_2           | (2,k2)
   *
   * the merged indices will be
   *
   * function in {f,g}| index
   * -----------------|------
   *  f_0             | (0,i0)
   *  f_1             | (0,i1)
   *  f_2             | (1,i2)
   *  g_0             | (2,k0)
   *  g_1             | (3,k1)
   *  g_2             | (4,k2)
   */
  struct FlatLexicographic
    : public IndexMergingStrategy
  {};

  /**
   * \brief Interleaved merging of direct children without blocking.
   *
   * \ingroup FunctionSpaceBasesUtilities
   *
   * Example: For two children {f} and {g} with multi-indices (all i* can be multi-indices, themselves)
   *
   * function in {f}| index
   * ---------------|------
   *  f_0           |(0,i0)
   *  f_1           |(1,i1)
   *  f_2           |(2,i2)
   *
   * function in {g}| index
   * ---------------|------
   *  g_0           | (0,i0)
   *  g_1           | (1,i1)
   *  g_2           | (2,i2)
   *
   * the merged indices will be
   *
   * function in {f,g}| index
   * -----------------|------
   *  f_0             | (0,i0)
   *  g_0             | (1,i0)
   *  f_1             | (2,i1)
   *  g_1             | (3,i1)
   *  f_2             | (4,i2)
   *  g_2             | (5,i2)
   */
  struct FlatInterleaved
    : public IndexMergingStrategy
  {};

  /**
   * \brief Lexicographic merging of direct children with blocking (i.e. creating one block per direct child).
   *
   * \ingroup FunctionSpaceBasesUtilities
   *
   * Example: For two children {f} and {g} with multi-indices (all i*,k* can be multi-indices, themselves)
   *
   * function in {f}| index
   * ---------------|------
   *  f_0           |(i0)
   *  f_1           |(i1)
   *  f_2           |(i2)
   *
   * function in {g}| index
   * ---------------|------
   *  g_0           | (k0)
   *  g_1           | (k1)
   *  g_2           | (k2)
   *
   * the merged indices will be
   *
   * function in {f,g}| index
   * -----------------|------
   *  f_0             | (0,i0)
   *  f_1             | (0,i1)
   *  f_2             | (0,i2)
   *  g_0             | (1,k0)
   *  g_1             | (1,k1)
   *  g_2             | (1,k2)
   */
  struct BlockedLexicographic
    : public IndexMergingStrategy
  {};

  /**
   * \brief Interleaved merging of direct children with blocking (i.e. creating blocks at the leaves containing one leaf per child each).
   *
   * \ingroup FunctionSpaceBasesUtilities
   *
   * Example: For two children {f} and {g} with multi-indices (all i* can be multi-indices, themselves)
   *
   * function in {f}| index
   * ---------------|------
   *  f_0           |(i0)
   *  f_1           |(i1)
   *  f_2           |(i2)
   *
   * function in {g}| index
   * ---------------|------
   *  g_0           | (i0)
   *  g_1           | (i1)
   *  g_2           | (i2)
   *
   * the merged indices will be
   *
   * function in {f,g}| index
   * -----------------|------
   *  f_0             | (i0,0)
   *  g_0             | (i0,1)
   *  f_1             | (i1,0)
   *  g_1             | (i1,1)
   *  f_2             | (i2,0)
   *  g_2             | (i2,1)
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
