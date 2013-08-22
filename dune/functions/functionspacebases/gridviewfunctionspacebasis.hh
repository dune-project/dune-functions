
template<typename T>
struct FunctionSpaceBasisTraits;
{
  typedef calc<T>::LocalBasisView;
};


template<typename GV, ...>
class GridViewFunctionSpaceBasis
{

  const GridView& gridView() const;

  /** maximum local size for any element on the GridView */
  size_type maxLocalSize() const;

};

//! TypeTree node
template<typename GridViewFunctionSpaceBasis>
class GridViewLocalBasisView
{

  //! Bind to element - only callable on root node
  void bind(const Element& e);

  //! Unbind from element - only hint, only callable on root node
  void unbind();

  some_type& tree() const;

};


template<typename GridViewFunctionSpaceBasis>
class GridViewLocalBasisViewTreeNode
{

  typedef typename GridViewFunctionSpaceBasis::Element Element;

  //! Return current element, throw if unbound
  const Element& element() const;

  //! size of subtree rooted in this node (element-local)
  size_type size() const;

  //! size of complete tree (element-local)
  size_type rootSize() const;

  //! Maps from subtree index set [0..size-1] into root index set (element-local) [0..??-1]
  size_type localIndex(size_type i) const;

  //! Maps from subtree index set [0..size-1] to a globally unique DOF index in global basis (pair of multi-indices)
  const DOFIndex& dofIndex(size_type i) const;

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type maxSize() const;

  //! maximum size of complete tree for any element of the global basis
  size_type maxRootSize() const;

  //! to be discussed
  const GridViewFunctionSpaceBasis& globalBasis() const;

};


template<typename GridViewFunctionSpaceBasis, typename... Children>
class CompositeGridViewLocalBasisViewTreeNode
  : public GridViewLocalBasisViewTreeNode<GridViewFunctionSpaceBasis>
  , public TypeTree::CompositeNode<Children...>
{};

template<typename GridViewFunctionSpaceBasis, typename Child, size_type n>
class PowerGridViewLocalBasisViewTreeNode
  : public GridViewLocalBasisViewTreeNode<GridViewFunctionSpaceBasis>
  , public TypeTree::PowerNode<Child,n>
{};


template<typename GridViewFunctionSpaceBasis>
class LeafGridViewLocalBasisViewTreeNode
  : public GridViewLocalBasisViewTreeNode<GridViewFunctionSpaceBasis>
  , public TypeTree::LeafNode
{

  FiniteElement finiteElement() const;

  //! Generate DOF indices for current subtree into range starting at it
  //! \param it iterator over a container of DOFIndex
  //! \return iterator past the last written element (STL-style)
  template<typename DOFIndexIterator>
  DOFIndexIterator generateDOFIndices(DOFIndexIterator it) const;

};
