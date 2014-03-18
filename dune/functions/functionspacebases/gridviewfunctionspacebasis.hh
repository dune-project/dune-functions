// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH

template<typename T>
struct FunctionSpaceBasisTraits;
{
  typedef calc<T>::LocalBasisView;
};


template<typename GV, ...>
class GridViewFunctionSpaceBasis
{
  typedef std::size_t size_type;
  typedef typename Dune::ReservedVector<size_type, 42> MultiIndex;

  /**
   * \brief
   */
  const GridView& gridView() const;

  /**
   * \brief maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * max{GridViewLocalBasisView(e).tree().size() | e in GridView}
   */
  size_type maxLocalSize() const;

  //! Return number possible values for next position in multi index
  size_type subIndexCount(const MultiIndex& index) const;

  /**
   * \brief Return local view for basis
   *
   * Perhaps we must move the construction outside
   * of the global basis in order to calm the compiler
   * when instanciating the TMP constructing the local view.
   */
  GridViewLocalBasisView localView() const;
};

//! TypeTree node
template<typename GridViewFunctionSpaceBasis>
class GridViewLocalBasisView
{
  typedef typename GridViewFunctionSpaceBasis::size_type size_type;
  typedef typename GridViewFunctionSpaceBasis::MultiIndex MultiIndex;

  //! Bind to element
  void bind(const Element& e);

  //! Unbind from element - only hint
  void unbind();

  /**
   * \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns some_type // This is tree
   */
  some_type& tree() const;
};


template<typename GridViewFunctionSpaceBasis>
class GridViewLocalBasisViewTreeNode
{
  typedef typename GridViewFunctionSpaceBasis::size_type size_type;
  typedef typename GridViewFunctionSpaceBasis::MultiIndex MultiIndex;

  typedef typename GridViewFunctionSpaceBasis::Element Element;

  //! Return current element, throw if unbound
  const Element& element() const;

  //! size of subtree rooted in this node (element-local)
//  size_type subTreeSize() const;

  //! maximum size of subtree rooted in this node for any element of the global basis
//  size_type maxSubTreeSize() const;

  //! size of complete tree (element-local)
  size_type localSize() const;

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const;

  //! maximum size of complete tree for any element of the global basis
  size_type maxLocalSize() const;

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis (pair of multi-indices)
  const MultiIndex& globalIndex(size_type i) const;

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

  //! Generate multi indices for current subtree into range starting at it
  //! \param it iterator over a container of MultiIndex
  //! \return iterator past the last written element (STL-style)
  template<typename MultiIndexIndexIterator>
  MultiIndexIndexIterator generateMultiIndexIndices(MultiIndexIndexIterator it) const;

};

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH
