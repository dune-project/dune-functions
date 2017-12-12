// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TREEDATA_HH
#define DUNE_FUNCTIONS_COMMON_TREEDATA_HH

#include <type_traits>
#include <utility>

#include <dune/typetree/traversal.hh>

namespace Dune {
namespace Functions {

/**
* \brief Container allowing to attach data to each node of a tree
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
class TreeData
{
public:
  /// Type of tree the data is associated with
  using Tree = T;

  /// Type used for indices and size information
  using size_type = typename Tree::size_type;

  /// Set if data should only be associated to the leafs
  static const bool leafOnly = LO;

  /// Template to determine the data type for given node type
  template<class Node>
  using NodeData = ND<Node>;

public:
  /// Default constructor
  TreeData()
    : tree_(nullptr)
  {}

  /**
    * \brief Construct from tree
    *
    * This default creates the data object associated to each node in the tree.
    * A reference to the tree is stored because it's needed for destruction
    * of the tree data.
    * See also \ref init.
    **/
  TreeData(const Tree& tree)
    : tree_(&tree)
  {
    initData();
  }

  /// Copy constructor
  TreeData(const TreeData& other)
    : TreeData(*other.tree_)
  {
    copyData(other);
  }

  /// Move constructor
  TreeData(TreeData&& other)
    : TreeData()
  {
    swap(other);
  }

  /**
    * \brief Initialize from tree
    *
    * This default creates the data object associated to each node in the tree.
    * A reference to the tree is stored because it's needed for destruction
    * of the tree data.
    **/
  void init(const Tree& tree)
  {
    destroyData();
    tree_ = &tree;
    initData();
  }

  /// Copy and Move assignment
  TreeData& operator=(TreeData other)
  {
    swap(other);
    return *this;
  }

  /// Destructor
  ~TreeData()
  {
    destroyData();
  }

  /// Return the attached tree (maybe nullptr)
  Tree const* tree() const
  {
    return tree_;
  }

  /// Get mutable reference to data associated to given node
  template <class Node>
  NodeData<Node>& operator[](const Node& node)
  {
    return *(NodeData<Node>*)(data_[node.treeIndex()]);
  }

  /// Get reference to data associated to given node
  template <class Node>
  const NodeData<Node>& operator[](const Node& node) const
  {
    return *(NodeData<Node>*)(data_[node.treeIndex()]);
  }

  /// Swap tree and data container with `other`
  void swap(TreeData& other)
  {
    using std::swap;
    swap(tree_, other.tree_);
    swap(data_, other.data_);
  }

protected:
  // For each node allocate memory and store the void* in the \ref data_
  void initData()
  {
    std::size_t s = 0;
    apply([&s](const auto& node, auto) {
      s = std::max(s, node.treeIndex()+1);
    });

    data_.resize(s, nullptr);
    apply([this](const auto& node, auto) {
      using Node = std::remove_reference_t<decltype(node)>;
      data_[node.treeIndex()] = new NodeData<Node>;
    });
  }

  // Deep copy of node data
  void copyData(const TreeData& other)
  {
    apply([&other,this](const auto& node, auto) {
      (*this)[node] = other[node];
    });
  }

  // For each node, delete the allocated node data
  void destroyData()
  {
    apply([this](const auto& node, auto) {
      using Node = std::remove_reference_t<decltype(node)>;
      auto* p = (NodeData<Node>*)(data_[node.treeIndex()]);
      delete p;
    });
    tree_ = nullptr;
  }

protected:
  // apply functor `func` to each (leaf) node of the tree, but only if `tree_` is set
  template <class Func>
  void apply(Func&& func)
  {
    if (tree_)
      applyImpl(func, std::integral_constant<bool, leafOnly>{});
  }

  template <class Func>
  void applyImpl(Func&& func, std::true_type) { forEachLeafNode(*tree_, func); }

  template <class Func>
  void applyImpl(Func&& func, std::false_type) { forEachNode(*tree_, func); }

protected:
  const Tree* tree_ = nullptr;
  std::vector<void*> data_;
};



} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_TREEDATA_HH
