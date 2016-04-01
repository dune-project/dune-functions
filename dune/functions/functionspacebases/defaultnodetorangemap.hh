// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTNODETORANGEMAP_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTNODETORANGEMAP_HH


#include <dune/common/concept.hh>

#include <dune/functions/functionspacebases/concepts.hh>

#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>


namespace Dune {
namespace Functions {



/**
 * \brief A simple node to range map using lexicographic ordering
 *
 * On construction this map is provided with an ansatz subtree
 * of the tree in an actual basis. The operator() then maps a
 * a pair (node,y) of a leaf node in this subtree and an object
 * y of type range to y[i] where i is the lexicographic index
 * of the node in the subtree.
 *
 * This allows to assoziate each leaf node in the subtree
 * with a component in the range type which is necessary
 * when interpolating a vector field into a nontrivial
 * subtree.
 *
 * Notice that the lexicographic index is only computed
 * wrt the leaf nodes.
 */
template<class Tree>
struct DefaultNodeToRangeMap
{

  // A simple visitor for computing lexicographic
  // subtree indices. To identify a leaf node
  // we use its treeIndex() which is unique
  // wrt the whole tree and store the computed
  // index in a vector indexed by the tree indices.
  struct Visitor
    : public TypeTree::TreeVisitor
    , public TypeTree::DynamicTraversal
  {
    Visitor(std::vector<std::size_t>& indices) :
      indices_(indices),
      counter_(0)
    {}

    template<typename Node, typename TreePath>
    void leaf(Node& node, TreePath treePath)
    {
      if (indices_.size() < node.treeIndex()+1)
        indices_.resize(node.treeIndex()+1);
      indices_[node.treeIndex()] = counter_;
      ++counter_;
    }

    std::vector<std::size_t>& indices_;
    std::size_t counter_;
  };

  /**
   * \brief Construct DefaultNodeToRangeMap
   *
   * \param tree The tree needed to initialize the map
   *
   * This only uses tree.treeIndex() such that the
   * tree can be in unbound state. Furthermore
   * no reference or pointer to the tree is stored
   * and passing a temporary is OK here.
   */
  DefaultNodeToRangeMap(const Tree& tree)
  {
    TypeTree::applyToTree(tree, Visitor(indices_));
  }

  template<class Node, class Range,
    typename std::enable_if<
      models<Concept::HasIndexAccess, Range, decltype(std::declval<Node>().treeIndex())>() and not Tree::isLeaf, int>::type = 0>
  auto operator()(const Node& node, Range&& y) const
    -> decltype(y[0])
  {
    return y[indices_[node.treeIndex()]];
  }

  template<class Node, class Range,
    typename std::enable_if< not models<Concept::HasIndexAccess, Range, decltype(std::declval<Node>().treeIndex())>() or Tree::isLeaf, int>::type = 0>
  auto operator()(const Node& node, Range&& y) const
    -> decltype(std::forward<Range>(y))
  {
    return std::forward<Range>(y);
  }

  std::vector<std::size_t> indices_;
  std::size_t counter_;
};



template<class Tree>
DefaultNodeToRangeMap<Tree> makeDefaultNodeToRangeMap(const Tree& tree)
{
  return DefaultNodeToRangeMap<Tree>(tree);
}



template<class Basis, class TreePath>
auto makeDefaultNodeToRangeMap(const Basis& basis, TreePath&& treePath)
  -> decltype(makeDefaultNodeToRangeMap(TypeTree::child(basis.localView().tree(),treePath)))
{
  auto&& localView = basis.localView();
  localView.bind(*basis.gridView().template begin<0>());
  auto&& tree = TypeTree::child(localView.tree(),treePath);
  return makeDefaultNodeToRangeMap(tree);
}



} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTNODETORANGEMAP_HH
