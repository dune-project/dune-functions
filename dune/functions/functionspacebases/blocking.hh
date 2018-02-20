// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BLOCKING_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BLOCKING_HH

#include <tuple>
#include <utility>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

namespace Dune { namespace Functions
{
  // Lightweight representation of (hierarchic) block structure
  namespace Blocking
  {
    struct Unknown {};

    /// \brief Blocking structure corresponding to IndexMergingStrategies \ref FlatLexicographic and \ref FlatInterleaved
    struct Flat
    {
      template <class Index>
      Flat operator[](const Index&) const { return {}; }
    };

    /// \brief Blocking structure corresponding to IndexMergingStrategy \ref LeafBlockedInterleaved
    template <std::size_t N>
    struct LeafBlocked
    {
      template <class Index>
      Flat operator[](const Index&) const { return {}; }
    };

    /// \brief Blocking structure corresponding to IndexMergingStrategy \ref BlockedLexicographic
    template <class Tag0, class... Tags>
    struct Blocked
    {
      template <std::size_t i>
      auto operator[](const index_constant<i>&) const { return std::get<i>(t); }

      Tag0 operator[](std::size_t /*i*/) const { return {}; }

      std::tuple<Tag0,Tags...> t;
    };

  } // end namespace tag


  namespace Concept
  {
    template <class T>
    struct IsFlat : std::false_type {};

    template <>
    struct IsFlat<Blocking::Flat> : std::true_type {};

    template <class T>
    struct IsBlocked : std::false_type {};

    template <class... Tags>
    struct IsBlocked<Blocking::Blocked<Tags...>> : std::true_type {};

    template <class T>
    static constexpr bool Flat() { return IsFlat<T>::value; }

    template <class T>
    static constexpr bool Blocked() { return IsBlocked<T>::value; }

  } // end namespace Concept


  namespace Impl
  {
    // Handle non-leaf nodes
    template <template <class...> class Builder, class NodeFactory>
    struct CompositeBlocking
    {
      using type = Blocking::Flat;
    };


    // Type is unknown if 1. level flat, 2. level blocked
    template <bool /*false*/, class... T>
    struct FlatTag { using type = Blocking::Unknown; };

    // Type is flat if all sub-Blockings are flat
    template <class... T>
    struct FlatTag<true, T...> { using type = Blocking::Flat; };

    /// \brief collapse flat tags to flat, if hierarchy is all flat, otherwise to tag::unknown
    template <class... T>
    using FlatTag_t = typename FlatTag<Std::conjunction<std::is_same<Blocking::Flat,T>...>::value, T...>::type;


    // Flat index-merging strategy
    template <class NodeFactory, class IndexMergingStrategy>
    struct BlockingImpl
    {
      using type = typename Impl::CompositeBlocking<FlatTag_t,NodeFactory>::type;
    };

    // Block index-merging strategy
    template <class NF>
    struct BlockingImpl<NF, BasisBuilder::BlockedLexicographic>
    {
      using type = typename Impl::CompositeBlocking<Blocking::Blocked, NF>::type;
    };


    template <class NF>
    struct Children
        : index_constant<0> {};

    template <class MI, class IMS, class... SF>
    struct Children<CompositeNodeFactory<MI, IMS, SF...>>
        : index_constant<sizeof...(SF)> {};

    template <class MI, class IMS, class SF, std::size_t C>
    struct Children<PowerNodeFactory<MI, IMS, SF, C>>
        : index_constant<C> {};


    // Leaf-Block index-merging strategy
    template <class NF>
    struct BlockingImpl<NF, BasisBuilder::LeafBlockedInterleaved>
    {
      using type = Blocking::LeafBlocked<Children<NF>::value>;
    };

  } // end namespace Impl


  // Handle leaf-nodes in basis-tree
  template <class NodeFactory, class = void>
  struct BlockingNodeFactory
  {
    using type = Blocking::Flat;
  };

  // Handle all nodes with an index-merging trategy
  template <class NF>
  struct BlockingNodeFactory<NF, void_t<typename NF::IndexMergingStrategy>>
  {
    using type = typename Impl::BlockingImpl<NF, typename NF::IndexMergingStrategy>::type;
  };

  /// \brief Extract block structure of NodeFactory
  template <class NodeFactory>
  using BlockingNodeFactory_t = typename BlockingNodeFactory<NodeFactory>::type;

  /// \brief Extract blocking structure of GlobalBasis
  template <class GlobalBasis>
  using Blocking_t = BlockingNodeFactory_t<typename GlobalBasis::NodeFactory>;


  namespace Impl
  {
    template <template <class...> class Builder, class MI, class IMS, class... SF>
    struct CompositeBlocking<Builder, CompositeNodeFactory<MI, IMS, SF...>>
    {
      using type = Builder<BlockingNodeFactory_t<SF>...>;
    };

    template <template <class...> class Builder, class MI, class IMS, class SF, std::size_t C>
    struct CompositeBlocking<Builder, PowerNodeFactory<MI, IMS, SF, C>>
    {
      template <std::size_t I> using SubFactories = SF;
      template <class Idx>        struct Expand;
      template <std::size_t... I> struct Expand<std::index_sequence<I...>>
      {
        using type = Builder<BlockingNodeFactory_t<SubFactories<I>>...>;
      };

      using type = typename Expand<std::make_index_sequence<C>>::type;
    };

  } // end namespace Impl

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BLOCKING_HH
