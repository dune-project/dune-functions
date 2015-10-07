// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH


#include <dune/functions/common/concept.hh>

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



} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
