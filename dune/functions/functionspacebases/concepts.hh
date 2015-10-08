// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH


#include <dune/functions/common/concept.hh>
#include <dune/functions/common/utility.hh>

#include <dune/functions/functionspacebases/nodes.hh>


namespace Dune {
namespace Functions {
namespace Concept {



struct HasResize
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.resize(0)
  );
};



struct HasSizeMethod
{
  template<class C>
  auto require(C&& c) -> decltype(
    c.size()
  );
};



struct HasIndexAccess
{
  template<class C, class I>
  auto require(C&& c, I&& i) -> decltype(
    c[i]
  );
};


// Concept for a BasisNode in a local ansatz tree
struct BasisNode
{
  template<class N>
  auto require(N&& node) -> decltype(
    requireType<typename N::size_type>(),
    requireType<typename N::TreePath>(),
    requireConvertible<typename N::size_type>(node.size()),
    requireConvertible<typename N::size_type>(node.offset()),
    requireConvertible<typename N::size_type>(node.localIndex(std::declval<typename N::size_type>())),
    requireConvertible<typename N::size_type>(node.treeIndex()),
    requireConvertible<typename N::TreePath>(node.treePath()),
    requireBaseOf<BasisNodeMixin<typename N::size_type, typename N::TreePath>, N>()
  );
};



// Concept for a LeafBasisNode in a local ansatz tree
template<class GridView>
struct LeafBasisNode : Refines<BasisNode>
{
  template<class N>
  auto require(N&& node) -> decltype(
    requireType<typename N::Element>(),
    requireType<typename N::FiniteElement>(),
    requireConvertible<typename N::Element>(node.element()),
    requireConvertible<const typename N::FiniteElement&>(node.finiteElement()),
    requireConvertible<typename N::Element>(*(std::declval<GridView>().template begin<0>())),
    requireBaseOf<Dune::Functions::LeafBasisNode<typename N::size_type, typename N::TreePath>, N>()
  );
};


template<class GridView>
struct BasisTree;

// Concept for a PowerBasisNode in a local ansatz tree
template<class GridView>
struct PowerBasisNode : Refines<BasisNode>
{
  template<class N>
  auto require(N&& node) -> decltype(
    requireBaseOf<Dune::Functions::PowerBasisNode<typename N::size_type, typename N::TreePath, typename N::ChildType, N::CHILDREN>, N>(),
    requireConcept<BasisTree<GridView>, typename N::ChildType>()
  );
};


// Concept for a CompositeBasisNode in a local ansatz tree
template<class GridView>
struct CompositeBasisNode : Refines<BasisNode>
{
  template<class ST, class TP>
  struct FixArgs
  {
    template<class...T>
    using CompositeBasisNode = typename Dune::Functions::CompositeBasisNode<ST, TP, T...>;
  };

  template<class N>
  auto require(const N& node) -> decltype(
    requireBaseOf<ExpandTuple<FixArgs<typename N::size_type, typename N::TreePath>::template CompositeBasisNode, typename N::ChildTypes>, N>(),
    requireConceptForTupleEntries<BasisTree<GridView>, typename N::ChildTypes>()
  );
};


// Concept for a full local BasisTree
template<class GridView>
struct BasisTree : Refines<BasisNode>
{
  template<class N>
  auto require(const N& node) -> decltype(
    requireConcept<typename std::conditional< N::isLeaf, LeafBasisNode<GridView>, BasisNode>::type, N>(),
    requireConcept<typename std::conditional< N::isPower, PowerBasisNode<GridView>, BasisNode>::type, N>(),
    requireConcept<typename std::conditional< N::isComposite, CompositeBasisNode<GridView>, BasisNode>::type, N>()
  );
};


// Concept for a NodeIndexSet
template<class NodeFactory>
struct NodeIndexSet
{
  template<class I>
  auto require(const I& indexSet) -> decltype(
    requireType<typename I::size_type>(),
    requireType<typename I::MultiIndex>(),
    requireType<typename I::NodeFactory>(),
    requireType<typename I::Node>(),
    requireSameType<typename I::NodeFactory, NodeFactory>(),
    const_cast<I&>(indexSet).bind(std::declval<typename I::Node>()),
    const_cast<I&>(indexSet).unbind(),
    requireConvertible<typename I::size_type>(indexSet.size()),
    requireConvertible<typename I::MultiIndex>(indexSet.index(std::declval<typename I::size_type>()))
  );
};


// Concept for a NodeFactory
template<class GridView>
struct NodeFactory
{
  using RootTreePath = decltype(TypeTree::hybridTreePath());

  template<class F>
  auto require(const F& factory) -> decltype(
    requireType<typename F::GridView>(),
    requireType<typename F::size_type>(),
    requireType<typename F::MultiIndex>(),
    requireType<typename F::SizePrefix>(),
    requireType<typename F::template Node<RootTreePath>>(),
    requireType<typename F::template IndexSet<RootTreePath>>(),
    requireSameType<typename F::GridView, GridView>(),
    const_cast<F&>(factory).initializeIndices(),
    requireConvertible<typename F::GridView>(factory.gridView()),
    requireConvertible<typename F::template Node<RootTreePath>>(factory.node(RootTreePath())),
    requireConvertible<typename F::template IndexSet<RootTreePath>>(factory.template indexSet<RootTreePath>()),
    requireConvertible<typename F::size_type>(factory.size()),
    requireConvertible<typename F::size_type>(factory.size(std::declval<typename F::SizePrefix>())),
    requireConvertible<typename F::size_type>(factory.dimension()),
    requireConvertible<typename F::size_type>(factory.maxNodeSize()),
    requireConcept<BasisTree<typename F::GridView>>(factory.node(RootTreePath())),
    requireConcept<NodeIndexSet<F>>(factory.template indexSet<RootTreePath>())
  );
};


} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
