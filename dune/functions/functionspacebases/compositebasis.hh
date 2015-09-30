// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH

#include <tuple>

#include <dune/common/reservedvector.hh>

#include <dune/typetree/compositenode.hh>
#include <dune/typetree/utility.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/staticforloop.hh>
#include <dune/functions/functionspacebases/basistags.hh>



namespace Dune {
namespace Functions {


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
 * This node factory represente a composition of several given node factories.
 * Its node type is a CompositeBasisNodes for the given subnodes.
 *
 * \tparam MI Type to be used for multi-indices
 * \tparam IT A tag describing how global indices are build
 * \tparam SF The sub-node factories
 */
template<class MI, class IT, class... SF>
class CompositeNodeFactory
{
  static const std::size_t children = sizeof...(SF);

  template<class, class, class, class...>
  friend class CompositeNodeIndexSet;

public:

  template<std::size_t k>
  using SubFactory = typename std::tuple_element<k, std::tuple<SF...>>::type;

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename SubFactory<0>::GridView;
  using size_type = typename SubFactory<0>::size_type;
  using IndexTag = IT;

  template<class TP, std::size_t k>
  using SubNode = typename SubFactory<k>::template Node<decltype(TypeTree::push_back(TP(), TypeTree::index_constant<k>()))>;

  template<class TP, std::size_t k>
  using SubIndexSet = typename SubFactory<k>::template IndexSet<decltype(TypeTree::push_back(TP(), TypeTree::index_constant<k>()))>;

  template<class, class>
  struct NodeHelper;

  template<class TP, std::size_t... k>
  struct NodeHelper<TP, TypeTree::Std::index_sequence<k...>>
  {
    using type = CompositeBasisNode<size_type, TP, SubNode<TP, k>... >;
  };

  template<class TP>
  using Node = typename NodeHelper<TP, TypeTree::Std::make_index_sequence<sizeof...(SF)>>::type;

  template<class TP>
  using IndexSet = CompositeNodeIndexSet<MI, TP, IT, SF...>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, MultiIndex::max_size()+1>;

private:

  using SubMultiIndex = MI;

public:

  /** \brief Constructor for a given grid view object */

  template<class... SFArgs>
  CompositeNodeFactory(SFArgs&&... sfArgs) :
    subFactories_(std::forward<SFArgs>(sfArgs)...)
  {}


  struct Lambda_initializeIndices
  {
    template<class I, class SFT>
    void operator()(const I& i, SFT& subFactories)
    {
      std::get<I::value>(subFactories).initializeIndices();
    }
  };

  void initializeIndices()
  {
    staticForLoop<0, sizeof...(SF)>(Lambda_initializeIndices(), subFactories_);
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return std::get<0>(subFactories_).gridView();
  }

  struct Lambda_node
  {
    template<class I, class Node, class SFT, class TP>
    void operator()(const I& i, Node& node, const SFT& subFactories, const TP& tp)
    {
      node.template setChild<I::value>(std::make_shared<SubNode<TP, I::value>>(std::get<I::value>(subFactories).node(TypeTree::push_back(tp, i))));
    }
  };

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    auto node = Node<TP>(tp);
    staticForLoop<0, sizeof...(SF)>(Lambda_node(), node, subFactories_, tp);
    return node;
  }

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

  struct Lambda_size
  {
    template<class I, class SFT>
    size_type operator()(const I&, const SFT& subFactories, const SizePrefix& prefix)
    {
      using SubFactory = typename std::tuple_element<I::value, SFT>::type;
      const SubFactory& subFactory = std::get<I::value>(subFactories);

      typename SubFactory::SizePrefix subPrefix;
      for(std::size_t i=1; i<prefix.size(); ++i)
        subPrefix.push_back(prefix[i]);
      return subFactory.size(subPrefix);
    }
  };

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return children;
    return forwardAsStaticIndex<children>(prefix[0], Lambda_size(), subFactories_, prefix);
  }

  struct Lambda_dimension
  {
    template<class I, class SFT>
    void operator()(const I& i, const SFT& subFactories, size_type& sum)
    {
      sum += std::get<i>(subFactories).dimension();
    }
  };

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    size_type r=0;
    staticForLoop<0, sizeof...(SF)>(Lambda_dimension(), subFactories_, r);
    return r;
  }

  struct Lambda_maxNodeSize
  {
    template<class I, class SFT>
    void operator()(const I& i, const SFT& subFactories, size_type& sum)
    {
      sum += std::get<i>(subFactories).maxNodeSize();
    }
  };

  size_type maxNodeSize() const
  {
    size_type r=0;
    staticForLoop<0, sizeof...(SF)>(Lambda_maxNodeSize(), subFactories_, r);
    return r;
  }

protected:
  std::tuple<SF...> subFactories_;
};



template<class MI, class TP, class IT, class... SF>
class CompositeNodeIndexSet
{
  static const std::size_t children = sizeof...(SF);

public:

  template<std::size_t k>
  using SubFactory = typename std::tuple_element<k, std::tuple<SF...>>::type;

  using GridView = typename SubFactory<0>::GridView;
  using size_type = typename SubFactory<0>::size_type;
  using IndexTag = IT;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = CompositeNodeFactory<MI, IT, SF...>;

  using Node = typename NodeFactory::template Node<TP>;

  template<std::size_t k>
  using SubTreePath = typename TypeTree::Child<Node,k>::TreePath;

  template<std::size_t k>
  using SubIndexSet = typename NodeFactory::template SubIndexSet<TP, k>;

  template<class>
  struct SubIndexSetTupleHelper;

  template<std::size_t... k>
  struct SubIndexSetTupleHelper<TypeTree::Std::index_sequence<k...>>
  {
    using type = std::tuple<SubIndexSet<k>...>;
  };

  using SubIndexSetTuple = typename SubIndexSetTupleHelper<TypeTree::Std::make_index_sequence<sizeof...(SF)>>::type;



  template<class SubFactoryTuple, size_t... I>
  static auto makeSubNodeIndexSetTuple(const SubFactoryTuple& subFactories, TypeTree::Std::index_sequence<I...>)
    -> decltype(std::make_tuple((std::get<I>(subFactories).template indexSet<SubTreePath<I>>())...))
  {
    return std::make_tuple((std::get<I>(subFactories).template indexSet<SubTreePath<I>>())...);
  }


  CompositeNodeIndexSet(const NodeFactory & nodeFactory) :
    nodeFactory_(&nodeFactory),
    subNodeIndexSetTuple_(makeSubNodeIndexSetTuple(nodeFactory_->subFactories_, TypeTree::Std::make_index_sequence<sizeof...(SF)>()))
  {}

  struct Lambda_bind
  {
    template<class I, class SNIT>
    void operator()(const I& i, SNIT& subNodeIndexSetTuple, const Node& node)
    {
      std::get<I::value>(subNodeIndexSetTuple).bind(node.template child<I::value>());
    }
  };

  void bind(const Node& node)
  {
    using namespace TypeTree::Indices;
    node_ = &node;
    staticForLoop<0, sizeof...(SF)>(Lambda_bind(), subNodeIndexSetTuple_, node);
  }

  struct Lambda_unbind
  {
    template<class I, class SNIT>
    void operator()(const I& i, SNIT& subNodeIndexSetTuple)
    {
      std::get<I::value>(subNodeIndexSetTuple).unbind();
    }
  };

  void unbind()
  {
    node_ = nullptr;
    staticForLoop<0, sizeof...(SF)>(Lambda_unbind(), subNodeIndexSetTuple_);
  }

  size_type size() const
  {
    return node_->size();
  }

  struct Lambda_index
  {
    template<class I, class SNIT>
    bool operator()(const I& i, SNIT& subNodeIndexSetTuple, size_type localIndex, size_type& offset, size_type& component, MultiIndex& multiIndex)
    {
      const auto& subNodeIndexSet = std::get<I::value>(subNodeIndexSetTuple);
      localIndex -= offset;
      size_type size = subNodeIndexSet.size();
      if (localIndex < size)
      {
        multiIndex = subNodeIndexSet.index(localIndex);
        component = i;
        return true;
      }
      offset += size;
      return false;
    }
  };

  MultiIndex index(size_type localIndex) const
  {
    size_type offset = 0;
    size_type component = 0;
    MultiIndex mi;
    staticFindInRange<0, sizeof...(SF)>(Lambda_index(), subNodeIndexSetTuple_, localIndex, offset, component, mi);
    mi.resize(mi.size()+1);

    for(std::size_t i=mi.size()-1; i>0; --i)
      mi[i] = mi[i-1];
    mi[0] = component;
    return mi;
  }

private:
  const NodeFactory* nodeFactory_;
  SubIndexSetTuple subNodeIndexSetTuple_;
  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<class IndexTag, class... SubFactoryTags>
struct CompositeNodeFactoryBuilder
{
  template<class MultiIndex, class GridView, class size_type=std::size_t>
  auto build(const GridView& gridView)
    -> CompositeNodeFactory<MultiIndex,  IndexTag, decltype(SubFactoryTags().template build<MultiIndex, GridView, size_type>(gridView))...>
  {
    return {SubFactoryTags().build<MultiIndex, GridView, size_type>(gridView)...};
  }
};

} // end namespace BasisBuilder::Imp

template<class IndexTag, class... SubFactoryTags>
Imp::CompositeNodeFactoryBuilder<IndexTag, SubFactoryTags...> composite(SubFactoryTags&&... tags, const IndexTag&)
{
  return{};
}

template<class... SubFactoryTags>
Imp::CompositeNodeFactoryBuilder<BasisTags::BlockFrontTag, SubFactoryTags...> composite(SubFactoryTags&&... tags)
{
  return{};
}

} // end namespace BasisBuilder



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_COMPOSITEBASIS_HH
