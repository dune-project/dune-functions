// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH


#include <dune/common/concept.hh>
#include <dune/common/reservedvector.hh>

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
    requireBaseOf<Dune::Functions::PowerBasisNode<typename N::ChildType, N::degree()>, N>(),
    requireConcept<BasisTree<GridView>, typename N::ChildType>()
  );
};

// Concept for a DynamicPowerBasisNode in a local ansatz tree
template<class GridView>
struct DynamicPowerBasisNode : Refines<BasisNode>
{
  template<class N>
  auto require(const N& node) -> decltype(
    requireBaseOf<Dune::Functions::DynamicPowerBasisNode<typename N::ChildType>, N>(),
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
  template<class N, class NodeTag = typename N::NodeTag>
  auto require(const N& node) -> decltype(
    requireConcept<std::conditional_t<N::isLeaf, LeafBasisNode<GridView>, BasisNode>, N>(),
    requireConcept<std::conditional_t<std::is_same_v<NodeTag, Dune::TypeTree::PowerNodeTag>,
                                      PowerBasisNode<GridView>, BasisNode>, N>(),
    requireConcept<std::conditional_t<std::is_same_v<NodeTag, Dune::TypeTree::DynamicPowerNodeTag>,
                                      DynamicPowerBasisNode<GridView>, BasisNode>, N>(),
    requireConcept<std::conditional_t<N::isComposite, CompositeBasisNode<GridView>, BasisNode>, N>()
  );
};


// Concept for a PreBasis
template<class GridView>
struct PreBasis
{
private:
  template<class PB>
  using MultiIndex = Dune::ReservedVector<typename PB::size_type, PB::multiIndexBufferSize>;

public:
  template<class PB>
  auto require(const PB& preBasis) -> decltype(
    requireType<typename PB::GridView>(),
    requireType<typename PB::size_type>(),
    requireType<typename PB::Node>(),
    requireConvertible<decltype(PB::maxMultiIndexSize), typename PB::size_type>(),
    requireConvertible<decltype(PB::maxMultiIndexSize), typename PB::size_type>(),
    requireConvertible<decltype(PB::multiIndexBufferSize), typename PB::size_type>(),
    requireTrue<PB::minMultiIndexSize <= PB::maxMultiIndexSize>(),
    requireTrue<PB::maxMultiIndexSize <= PB::multiIndexBufferSize>(),
    requireSameType<typename PB::GridView, GridView>(),
    const_cast<PB&>(preBasis).initializeIndices(),
    requireConvertible<typename PB::GridView>(preBasis.gridView()),
    requireConvertible<typename PB::Node>(preBasis.makeNode()),
    requireConvertible<typename PB::size_type>(preBasis.size()),
    requireConvertible<typename PB::size_type>(preBasis.size(std::declval<MultiIndex<PB>>())),
    requireConvertible<typename PB::size_type>(preBasis.dimension()),
    requireConvertible<typename PB::size_type>(preBasis.maxNodeSize()),
    requireSameType<decltype(const_cast<PB&>(preBasis).update(preBasis.gridView())),void>(),
    requireConcept<BasisTree<typename PB::GridView>>(preBasis.makeNode()),
    requireConvertible<typename std::vector<MultiIndex<PB>>::iterator>(
      preBasis.indices(
        preBasis.makeNode(),
        std::declval<typename std::vector<MultiIndex<PB>>::iterator>()))
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
    requireConvertible<bool>(localView.bound()),
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
