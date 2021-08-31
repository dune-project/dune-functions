// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TREEDATA_HH
#define DUNE_FUNCTIONS_COMMON_TREEDATA_HH


#warning This file is deprecated. Please use TreeContainer from dune-typetree instead.

#include <memory>

#include <dune/common/shared_ptr.hh>

#include <dune/typetree/pairtraversal.hh>

#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>

namespace Dune {
namespace Functions {

/**
 * \brief Mixin for visitors that should apply the same action on all nodes
 *
 * \deprecated
 *
 * \ingroup Utility
 *
 * By deriving from this you only have to implement apply(node,treepath)
 * in the derived class. This will be used for pre(...) and leaf(...).
 *
 * \tparam Derived Type of derived class implementing apply(node,treepath)
 * \tparam leafOnly Flag to enable leaf only traversal
 */
template<class SimpleNodeVisitorImp, bool leafOnly>
struct
[[deprecated("This is an implementation detail of the deprecated class TreeDate and thus deprecated itself.")]]
UniformNodeVisitor :
    public TypeTree::TreeVisitor,
    public TypeTree::DynamicTraversal
{
    // This is only enabled, if we want to incorporate inner nodes.
    // Checking leafOnly would be sufficient, but for SFINAE the
    // the enable_if condition must depend on the template parameter.
    template<typename Node, typename TreePath,
        typename std::enable_if<(not leafOnly) and (not Node::isLeaf), int>::type = 0>
    void pre(Node& node, TreePath treePath)
    {
      static_cast<SimpleNodeVisitorImp*>(this)->apply(node, treePath);
    }

    template<typename Node, typename TreePath,
        typename std::enable_if<(leafOnly) and (not Node::isLeaf), int>::type = 0>
    void pre(Node& node, TreePath treePath)
    {}

    template<typename Node, typename TreePath>
    void leaf(Node& node, TreePath treePath)
    {
      static_cast<SimpleNodeVisitorImp*>(this)->apply(node, treePath);
    }
};



/**
 * \brief Container allowing to attach data to each node of a tree
 *
 * \deprecated
 *
 * \ingroup Utility
 *
 * This provides operator[](Node) for accessing the data attached to the node.
 * For storing the data each node is identified via its treeIndex() method
 * which is supposed to return an index which is unique wrt the tree. These
 * indices need not to be consecutive but they are used to access an
 * internal vector<void*>. This may lead to wasted memory if the maximal
 * treeIndex() is much larger then the number of nodes within the tree.
 *
 * Before using the container it must be initialized by providing the
 * tree. The stored data objects will be created on initialization. Hence
 * the type of these data objects must be default constructible.
 *
 * Notice that the data per node can only be interpreted if the
 * node type is known. Hence the tree will be traversed on initilization,
 * copy, assignment, and destruction of a TreeData container.
 *
 * \tparam T Type of the tree
 * \tparam ND The data stored for a node of type Node will be of type ND<Node>
 * \tparam LO Set this flag if data should only be attached to leaf nodes.
 */
template<class T, template<class> class ND, bool LO>
class
[[deprecated("This class is deprecated. Please use TreeContainer from dune-typetree instead.")]]
TreeData
{

public:

  //! Type of tree the data is associated with
  using Tree = T;

  //! Type used for indices and size information
  using size_type = typename Tree::size_type;

  //! Set if data should only be associated to the leafs
  static const bool leafOnly = LO;

  //! Template to determine the data type for given node type
  template<class Node>
  using NodeData = ND<Node>;

protected:
  using RawContainer = std::vector<void*>;


  // Since we can generate the node data type only if
  // we know the type of the node, we have to do
  // initialization, copy, and destruction via a
  // tree traversal. Once we can use C++14 this can
  // be written in a much easier and more selfcontained
  // ways using generic lambda functions.
  // Until then we need explicite visitor classes for
  // each operation.

  struct InitVisitor :
    public UniformNodeVisitor<InitVisitor, leafOnly>
  {
    InitVisitor(RawContainer& data) :
      data_(data)
    {}

    template<typename Node, typename TreePath>
    void apply(Node& node, TreePath treePath)
    {
      auto&& index = node.treeIndex();
      if (data_.size() < index+1)
        data_.resize(index+1, nullptr);
      data_[index] = new NodeData<Node>;
    }


    RawContainer& data_;
  };

  struct DestroyVisitor :
    public UniformNodeVisitor<DestroyVisitor, leafOnly>
  {
    DestroyVisitor(RawContainer& data) :
      data_(data)
    {}

    template<typename Node, typename TreePath>
    void apply(Node& node, TreePath treePath)
    {
      auto&& index = node.treeIndex();
      auto p = (NodeData<Node>*)(data_[index]);
      delete p;
      data_[index] = nullptr;
    }

    RawContainer& data_;
  };

  struct CopyVisitor :
    public UniformNodeVisitor<CopyVisitor, leafOnly>
  {
    CopyVisitor(TreeData& thisTD, const TreeData& otherTD) :
      thisTD_(thisTD),
      otherTD_(otherTD)
    {}

    template<typename Node, typename TreePath>
    void apply(Node& node, TreePath treePath)
    {
      thisTD_[node] = otherTD_[node];
    }

    TreeData& thisTD_;
    const TreeData& otherTD_;
  };

public:

  //! Default constructor
  TreeData() :
    tree_(nullptr)
  {}

  /**
   * \brief Initialize from tree
   *
   * This default creates the data object associated to each node in the tree.
   * A reference to the tree is stored because it's needed for destruction
   * of the tree data.
   */
  void init(const Tree& tree)
  {
    if (tree_)
      destroy();
    tree_ = &tree;
    TypeTree::applyToTree(*tree_, InitVisitor(data_));
  }

  //! Copy constructor
  TreeData(const TreeData& other) :
    tree_(other.tree_)
  {
    TypeTree::applyToTree(*tree_, InitVisitor(data_));
    TypeTree::applyToTree(*tree_, CopyVisitor(*this, other));
  }

  //! Copy assignment
  TreeData& operator=(const TreeData& other)
  {
    if (tree_)
      TypeTree::applyToTree(*tree_, DestroyVisitor(data_));
    tree_ = other.tree_;
    TypeTree::applyToTree(*tree_, CopyVisitor(*this, other));
    return *this;
  }

  //! Destroy data
  void destroy()
  {
    if (tree_)
      TypeTree::applyToTree(*tree_, DestroyVisitor(data_));
    tree_ = nullptr;
  }

  //! Destructor
  ~TreeData()
  {
    if (tree_)
      TypeTree::applyToTree(*tree_, DestroyVisitor(data_));
  }

  //! Get mutable reference to data associated to given node
  template<class Node>
  NodeData<Node>& operator[](const Node& node)
  {
    return *(NodeData<Node>*)(data_[node.treeIndex()]);
  }

  //! Get reference to data associated to given node
  template<class Node>
  const NodeData<Node>& operator[](const Node& node) const
  {
    return *(NodeData<Node>*)(data_[node.treeIndex()]);
  }

protected:

  const Tree* tree_;
  RawContainer data_;
};



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_TREEDATA_HH
