// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH


#include <dune/common/concept.hh>

#include <dune/functions/common/utility.hh>

#include <dune/functions/functionspacebases/nodes.hh>


namespace Dune {
namespace Functions {
namespace Concept {

using namespace Dune::Concept;


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
  auto require(const N& node) -> decltype(
    requireType<typename N::size_type>(),
    requireConvertible<typename N::size_type>(node.size()),
    requireConvertible<typename N::size_type>(node.localIndex(std::declval<typename N::size_type>())),
    requireConvertible<typename N::size_type>(node.treeIndex()),
    requireBaseOf<BasisNodeMixin, N>()
  );
};



// Concept for a LeafBasisNode in a local ansatz tree
template<class GridView>
struct LeafBasisNode : Refines<BasisNode>
{
  template<class N>
  auto require(const N& node) -> decltype(
    requireType<typename N::Element>(),
    requireType<typename N::FiniteElement>(),
    requireConvertible<typename N::Element>(node.element()),
    requireConvertible<const typename N::FiniteElement&>(node.finiteElement()),
    requireSameType<typename N::Element, typename GridView::template Codim<0>::Entity>(),
    requireBaseOf<Dune::Functions::LeafBasisNode, N>()
  );
};


template<class GridView>
struct BasisTree;

// Concept for a PowerBasisNode in a local ansatz tree
template<class GridView>
struct PowerBasisNode : Refines<BasisNode>
{
  template<class N>
  auto require(const N& node) -> decltype(
    requireBaseOf<Dune::Functions::PowerBasisNode<typename N::ChildType, N::CHILDREN>, N>(),
    requireConcept<BasisTree<GridView>, typename N::ChildType>()
  );
};


// Concept for a CompositeBasisNode in a local ansatz tree
template<class GridView>
struct CompositeBasisNode : Refines<BasisNode>
{
  template<class N>
  auto require(const N& node) -> decltype(
    requireBaseOf<ExpandTuple<Dune::Functions::template CompositeBasisNode, typename N::ChildTypes>, N>(),
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
template<class PreBasis>
struct NodeIndexSet
{
  template<class I>
  auto require(const I& indexSet) -> decltype(
    requireType<typename I::size_type>(),
    requireType<typename I::MultiIndex>(),
    requireType<typename I::PreBasis>(),
    requireType<typename I::Node>(),
    requireSameType<typename I::PreBasis, PreBasis>(),
    const_cast<I&>(indexSet).bind(std::declval<typename I::Node>()),
    const_cast<I&>(indexSet).unbind(),
    requireConvertible<typename I::size_type>(indexSet.size()),
    requireConvertible<typename std::vector<typename I::MultiIndex>::iterator>(
      indexSet.indices(std::declval<typename std::vector<typename I::MultiIndex>::iterator>()))
  );
};


// Concept for a PreBasis
template<class GridView>
struct PreBasis
{
  template<class PB>
  auto require(const PB& preBasis) -> decltype(
    requireType<typename PB::GridView>(),
    requireType<typename PB::size_type>(),
    requireType<typename PB::MultiIndex>(),
    requireType<typename PB::SizePrefix>(),
    requireType<typename PB::Node>(),
    requireType<typename PB::IndexSet>(),
    requireSameType<typename PB::GridView, GridView>(),
    const_cast<PB&>(preBasis).initializeIndices(),
    requireConvertible<typename PB::GridView>(preBasis.gridView()),
    requireConvertible<typename PB::Node>(preBasis.makeNode()),
    requireConvertible<typename PB::IndexSet>(preBasis.makeIndexSet()),
    requireConvertible<typename PB::size_type>(preBasis.size()),
    requireConvertible<typename PB::size_type>(preBasis.size(std::declval<typename PB::SizePrefix>())),
    requireConvertible<typename PB::size_type>(preBasis.dimension()),
    requireConvertible<typename PB::size_type>(preBasis.maxNodeSize()),
    requireSameType<decltype(const_cast<PB&>(preBasis).update(preBasis.gridView())),void>(),
    requireConcept<BasisTree<typename PB::GridView>>(preBasis.makeNode()),
    requireConcept<NodeIndexSet<PB>>(preBasis.makeIndexSet())
  );
};



// Concept for a LocalView
template<class GlobalBasis>
struct LocalView
{
  template<class V>
  auto require(const V& localView) -> decltype(
    requireType<typename V::size_type>(),
    requireType<typename V::MultiIndex>(),
    requireType<typename V::GlobalBasis>(),
    requireType<typename V::Tree>(),
    requireType<typename V::GridView>(),
    requireType<typename V::Element>(),
    requireSameType<typename V::GlobalBasis, GlobalBasis>(),
    requireSameType<typename V::GridView, typename GlobalBasis::GridView>(),
    requireSameType<typename V::size_type, typename GlobalBasis::size_type>(),
    requireSameType<typename V::Element, typename GlobalBasis::GridView::template Codim<0>::Entity>(),
    const_cast<V&>(localView).bind(std::declval<typename V::Element>()),
    const_cast<V&>(localView).unbind(),
    requireConvertible<typename V::Tree>(localView.tree()),
    requireConvertible<typename V::size_type>(localView.size()),
    requireConvertible<typename V::MultiIndex>(localView.index(std::declval<typename V::size_type>())),
    requireConvertible<typename V::size_type>(localView.maxSize()),
    requireConvertible<typename V::GlobalBasis>(localView.globalBasis()),
    requireConcept<BasisTree<typename V::GridView>>(localView.tree()),
    0
  );
};



// Concept for a GlobalBasis
template<class GridView>
struct GlobalBasis
{
  template<class B>
  auto require(const B& basis) -> decltype(
    requireType<typename B::GridView>(),
    requireType<typename B::size_type>(),
    requireType<typename B::MultiIndex>(),
    requireType<typename B::SizePrefix>(),
    requireType<typename B::LocalView>(),
    requireSameType<typename B::GridView, GridView>(),
    requireConvertible<typename B::GridView>(basis.gridView()),
    requireConvertible<typename B::LocalView>(basis.localView()),
    requireConvertible<typename B::size_type>(basis.size()),
    requireConvertible<typename B::size_type>(basis.size(std::declval<typename B::SizePrefix>())),
    requireConvertible<typename B::size_type>(basis.dimension()),
    requireSameType<decltype(const_cast<B&>(basis).update(basis.gridView())),void>(),
    requireConcept<LocalView<B>>(basis.localView())
  );
};



} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
