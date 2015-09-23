// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH

#include <dune/typetree/leafnode.hh>
#include <dune/functions/functionspacebases/nodes.hh>

namespace Dune {
namespace Functions {



template<typename E, typename FE, typename ST, typename TP>
class GridFunctionSpaceBasisLeafNodeInterface :
  public LeafBasisNode<ST, TP>
{
  using Base = LeafBasisNode<ST, TP>;
public:
  typedef ST size_type;
  typedef E Element;
  typedef FE FiniteElement;

  using TreePath = typename Base::TreePath;

  GridFunctionSpaceBasisLeafNodeInterface(TreePath treePath = TreePath()) :
    Base(treePath)
  {}


  //! Return current element, throw if unbound
  virtual const Element& element() const = 0;

  virtual const FiniteElement& finiteElement() const = 0;

  //! size of subtree rooted in this node (element-local)
  virtual size_type size() const = 0;

};

template<typename GV, typename LV, typename IS, typename MI>
class GridViewFunctionSpaceBasis
{
public:

  typedef GV GridView;
  typedef std::size_t size_type;
  typedef LV LocalView;
  typedef IS IndexSet;
  typedef MI MultiIndex;

  /** \brief Obtain the grid view that the basis is defined on
   */
  virtual const GridView& gridView() const = 0;

  /**
   * \brief Return local view for basis
   *
   * Perhaps we must move the construction outside
   * of the global basis in order to calm the compiler
   * when instantiating the TMP constructing the local view.
   */
  virtual LocalView localView() const = 0;

};

#if 0
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

#endif

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GRIDVIEWFUNCTIONSPACEBASIS_HH
