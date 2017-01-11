// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH

#include <tuple>
#include <utility>

#include <dune/common/std/utility.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

#include <dune/functions/common/staticforloop.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>



namespace Dune {
namespace Functions {

namespace Imp {

  template<typename... T>
  struct SizeOf
    : public std::integral_constant<std::size_t,sizeof...(T)>
  {};

  template<typename... T>
  using index_sequence_for = std::make_index_sequence<SizeOf<T...>{}>;
}

// *****************************************************************************
// This is the reusable part of the composite bases. It contains
//
//   CompositeNodeFactory
//   CompositeNodeIndexSet
//
// The factory allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************


template<class MI, class TP, class IT, class... SF>
class CompositeNodeIndexSet;

/**
 * \brief A factory for composite bases
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This node factory represente a composition of several given node factories.
 * Its node type is a CompositeBasisNodes for the given subnodes.
 *
 * \tparam MI  Type to be used for global multi-indices
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child factories
 * \tparam SF  The child factories
 */
template<class MI, class IMS, class... SF>
class CompositeNodeFactory
{
public:

  //! Tuple of child factories
  using SubFactories = std::tuple<SF...>;

  //! The grid view that the FE basis is defined on
  using GridView = typename std::tuple_element<0, SubFactories>::type::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Strategy used to merge the global indices of the child factories
  using IndexMergingStrategy = IMS;

protected:
  static const std::size_t children = sizeof...(SF);

  template<class, class, class, class...>
  friend class CompositeNodeIndexSet;

  using ChildIndexTuple = IntegerSequenceTuple<Imp::index_sequence_for<SF...>>;

  template<class TP>
  struct FixedTP
  {

    template<class I>
    using IndexToSubTreePath = decltype(TypeTree::push_back(TP(), I()));

    using SubTreePaths = TransformTuple<IndexToSubTreePath, ChildIndexTuple>;

    template<class F, class SubTP>
    using FactoryToSubNode = typename F::template Node<SubTP>;

    using SubNodes = TransformTuple<FactoryToSubNode, SubFactories, SubTreePaths>;

    template<class F, class SubTP>
    using FactoryToSubIndexSet = typename F::template IndexSet<SubTP>;

    using SubIndexSets = TransformTuple<FactoryToSubIndexSet, SubFactories, SubTreePaths>;

    template<class... N>
    using SubNodesToNode = CompositeBasisNode<size_type, TP, N... >;

    using Node = ExpandTuple<SubNodesToNode, SubNodes>;
  };


public:

  //! Template mapping index of child to its factory type
  template<std::size_t k>
  using SubFactory = typename std::tuple_element<k, std::tuple<SF...>>::type;

  //! Template mapping root tree path to type of created tree node
  template<class TP>
  using Node = typename FixedTP<TP>::Node;

  //! Template mapping root tree path to type of created tree node index set
  template<class TP>
  using IndexSet = CompositeNodeIndexSet<MI, TP, IMS, SF...>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, MultiIndex::max_size()+1>;

  /**
   * \brief Constructor for given child factory objects
   *
   * The child factories will be stored as copies
   */
  template<class... SFArgs,
    disableCopyMove<CompositeNodeFactory, SFArgs...> = 0,
    enableIfConstructible<std::tuple<SF...>, SFArgs...> = 0>
  CompositeNodeFactory(SFArgs&&... sfArgs) :
    subFactories_(std::forward<SFArgs>(sfArgs)...)
  {
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    using namespace Dune::Hybrid;
    forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
      elementAt(subFactories_, i).initializeIndices();
    });
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return std::get<0>(subFactories_).gridView();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    using namespace Dune::Hybrid;
    forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
      elementAt(subFactories_, i).update(gv);
    });
  }

  /**
   * \brief Create tree node with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * By passing a non-trivial root tree path this can be used
   * to create a node suitable for being placed in a tree at
   * the position specified by the root tree path.
   */
  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    auto node = Node<TP>(tp);
    using namespace Dune::Hybrid;
    forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
      node.setChild( elementAt(subFactories_, i).node(TypeTree::push_back(tp, i)), i);
    });
    return node;
  }

  /**
   * \brief Create tree node index set with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * Create an index set suitable for the tree node obtained
   * by node(tp).
   */
  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return size({});
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return size(prefix, IndexMergingStrategy{});
  }

private:

  size_type size(const SizePrefix& prefix, BasisBuilder::BlockedLexicographic) const
  {
    if (prefix.size() == 0)
      return children;

    return Hybrid::switchCases(std::make_index_sequence<children>(), prefix[0], [&] (auto i) {
      const auto& subFactory = std::get<i.value>(subFactories_);
      typename std::decay<decltype(subFactory)>::type::SizePrefix subPrefix;
      for(std::size_t i=1; i<prefix.size(); ++i)
        subPrefix.push_back(prefix[i]);
      return subFactory.size(subPrefix);
    }, []() {
      return size_type(0);
    });
  }

  struct Lambda_size_flat_sizeInSubtree
  {
    template<class I, class SFT>
    size_type operator()(const I&, const SFT& subFactories, const SizePrefix& prefix, size_type& shiftedFirst, size_type& r)
    {
      using SubFactory = typename std::tuple_element<I::value, SFT>::type;
      const SubFactory& subFactory = std::get<I::value>(subFactories);
      if (shiftedFirst < subFactory.size())
      {
        typename SubFactory::SizePrefix subPrefix;
        subPrefix.push_back(shiftedFirst);
        for(std::size_t i=1; i<prefix.size(); ++i)
          subPrefix.push_back(prefix[i]);
        r = subFactory.size(subPrefix);
        return true;
      }
      shiftedFirst -= subFactory.size();
      return false;
    }
  };

  size_type size(const SizePrefix& prefix, BasisBuilder::FlatLexicographic) const
  {
    size_type r = 0;
    using namespace Dune::Hybrid;
    if (prefix.size() == 0)
      forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
        r += elementAt(subFactories_, i).size();
      });
    else {
      size_type shiftedFirst = prefix[0];
      staticFindInRange<0, sizeof...(SF)>(Lambda_size_flat_sizeInSubtree(), subFactories_, prefix, shiftedFirst, r);
    }
    return r;
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    size_type r=0;
    // Accumulate dimension() for all subfactories
    using namespace Dune::Hybrid;
    forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
      r += elementAt(subFactories_, i).dimension();
    });
    return r;
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    size_type r=0;
    // Accumulate maxNodeSize() for all subfactories
    using namespace Dune::Hybrid;
    forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
      r += elementAt(subFactories_, i).maxNodeSize();
    });
    return r;
  }

protected:
  std::tuple<SF...> subFactories_;
};



template<class MI, class TP, class IMS, class... SF>
class CompositeNodeIndexSet
{
  static const std::size_t children = sizeof...(SF);

public:

  template<std::size_t k>
  using SubFactory = typename std::tuple_element<k, std::tuple<SF...>>::type;

  using GridView = typename SubFactory<0>::GridView;
  using size_type = std::size_t;
  using IndexMergingStrategy = IMS;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = CompositeNodeFactory<MI, IMS, SF...>;

  using Node = typename NodeFactory::template Node<TP>;

  using SubTreePaths = typename NodeFactory::template FixedTP<TP>::SubTreePaths;
  using SubIndexSets = typename NodeFactory::template FixedTP<TP>::SubIndexSets;


  struct Lambda_FactoryToSubIndexSet
  {
    // transform a single (factory,subTreePath) pair to subIndexSet
    template<class Factory, class SubTP>
    auto operator()(const Factory& factory, const SubTP& subTP)
      ->decltype(factory.template indexSet<SubTP>())
    {
      return factory.template indexSet<SubTP>();
    }
  };

  CompositeNodeIndexSet(const NodeFactory & nodeFactory) :
    nodeFactory_(&nodeFactory),
    subNodeIndexSetTuple_(transformTuple(Lambda_FactoryToSubIndexSet(), nodeFactory_->subFactories_, SubTreePaths()))
  {}

  void bind(const Node& node)
  {
    node_ = &node;
    using namespace Dune::Hybrid;
    forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
      elementAt(subNodeIndexSetTuple_, i).bind(node.child(i));
    });
  }

  void unbind()
  {
    node_ = nullptr;
    using namespace Dune::Hybrid;
    forEach(Dune::Std::make_index_sequence<children>(), [&](auto i) {
      elementAt(subNodeIndexSetTuple_, i).unbind();
    });
  }

  size_type size() const
  {
    return node_->size();
  }

  MultiIndex index(size_type localIndex) const
  {
    return index(localIndex, IndexMergingStrategy{});
  }

  struct Lambda_index_flat
  {
    template<class I, class SNIT, class SFT>
    bool operator()(const I& i, SNIT& subNodeIndexSetTuple, const SFT& subFactories, size_type localIndex, size_type& rootOffset, MultiIndex& multiIndex)
    {
      const auto& subNodeIndexSet = std::get<I::value>(subNodeIndexSetTuple);
      size_type size = subNodeIndexSet.size();
      if (localIndex < size)
      {
        multiIndex = subNodeIndexSet.index(localIndex);
        multiIndex[0] += rootOffset;
        return true;
      }
      localIndex -= size;
      rootOffset += std::get<I::value>(subFactories).size();
      return false;
    }
  };

  MultiIndex index(const size_type& localIndex, BasisBuilder::FlatLexicographic) const
  {
    size_type shiftedLocalIndex = localIndex;
    size_type rootOffset = 0;
    MultiIndex mi;
    staticFindInRange<0, sizeof...(SF)>(Lambda_index_flat(), subNodeIndexSetTuple_, nodeFactory_->subFactories_, shiftedLocalIndex, rootOffset, mi);
    return mi;
  }

  struct Lambda_index
  {
    template<class I, class SNIT>
    bool operator()(const I& i, SNIT& subNodeIndexSetTuple, size_type& localIndex, size_type& component, MultiIndex& multiIndex)
    {
      const auto& subNodeIndexSet = std::get<I::value>(subNodeIndexSetTuple);
      size_type size = subNodeIndexSet.size();
      if (localIndex < size)
      {
        multiIndex = subNodeIndexSet.index(localIndex);
        component = i;
        return true;
      }
      localIndex -= size;
      return false;
    }
  };

  MultiIndex index(const size_type& localIndex, BasisBuilder::BlockedLexicographic) const
  {
    size_type shiftedLocalIndex = localIndex;
    size_type component = 0;
    MultiIndex mi;
    staticFindInRange<0, sizeof...(SF)>(Lambda_index(), subNodeIndexSetTuple_, shiftedLocalIndex, component, mi);
    mi.resize(mi.size()+1);

    for(std::size_t i=mi.size()-1; i>0; --i)
      mi[i] = mi[i-1];
    mi[0] = component;
    return mi;
  }

private:
  const NodeFactory* nodeFactory_;
  SubIndexSets subNodeIndexSetTuple_;
  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<class ST0>
constexpr std::size_t maxHelper(ST0&& i0)
{
  return i0;
}

template<class ST0, class... ST>
constexpr std::size_t maxHelper(ST0&& i0, ST&&... i)
{
  return (i0 > maxHelper(i...)) ? i0 : maxHelper(i...);
}

template<class IndexTag, class... SubFactoryTags>
struct CompositeNodeFactoryBuilder
{
  static const bool isBlocked = std::is_same<IndexTag,BlockedLexicographic>::value or std::is_same<IndexTag,LeafBlockedInterleaved>::value;

  static const std::size_t requiredMultiIndexSize=maxHelper(SubFactoryTags::requiredMultiIndexSize...) + (std::size_t)(isBlocked);

  template<class MultiIndex, class GridView>
  auto build(const GridView& gridView)
    -> CompositeNodeFactory<MultiIndex,  IndexTag, decltype(SubFactoryTags().template build<MultiIndex, GridView>(gridView))...>
  {
    return {SubFactoryTags().template build<MultiIndex, GridView>(gridView)...};
  }
};

template<class... Args>
auto compositeImp(std::tuple<Args...>)
  -> Imp::CompositeNodeFactoryBuilder<Args...>
{
  return {};
}

} // end namespace BasisBuilder::Imp



/**
 * \brief Create a factory builder that can build a CompositeNodeFactory
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam Args Types of child factory builders and IndexMergingStrategy type
 * \param args Child factory builder objects and an IndexMergingStrategy
 *
 * This is the overload used if the last argument is an IndexMergingStrategy.
 */
template<
  typename... Args,
  std::enable_if_t<Concept::isIndexMergingStrategy<typename LastType<Args...>::type>(),int> = 0>
auto composite(Args&&... args)
{
  return Imp::compositeImp(typename RotateTuple<Args...>::type{});
}

/**
 * \brief Create a factory builder that can build a CompositeNodeFactory
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam Args Types of child factory builders
 * \param args Child factory builder objects
 *
 * This is the overload used if no IndexMergingStrategy is supplied.
 * In this case the BasisBuilder::BlockedLexicographic strategy is used.
 */
template<
  typename... Args,
  std::enable_if_t<not Concept::isIndexMergingStrategy<typename LastType<Args...>::type>(),int> = 0>
auto composite(Args&&... args)
{
  return Imp::compositeImp(std::tuple<BasisBuilder::BlockedLexicographic,Args...>{});
}

} // end namespace BasisBuilder



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH
