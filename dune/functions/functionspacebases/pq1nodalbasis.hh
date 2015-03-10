// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH

#include <array>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
#include <dune/common/std/final.hh>
#else
 #ifndef DUNE_FINAL
  #define DUNE_FINAL
 #endif
#endif
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {


template<typename GV>
class PQ1NodalBasisLocalView;

template<typename GV>
class PQ1NodalBasisLeafNode;

/** \brief Nodal basis of a scalar first-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV>
class PQ1NodalBasis
: public GridViewFunctionSpaceBasis<GV,
                                    PQ1NodalBasisLocalView<GV>,
                                    std::array<std::size_t, 1> >
{
  static const int dim = GV::dimension;

public:
  typedef GV GridView;
  typedef std::size_t size_type;
  typedef PQ1NodalBasisLocalView<GV> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  /** \brief Constructor for a given grid view object */
  PQ1NodalBasis(const GridView& gv) :
    gridView_(gv)
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  /**
   * \brief maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * max{GridViewLocalBasisView(e).tree().size() | e in GridView}
   */
  size_type maxLocalSize() const DUNE_FINAL
  {
    return 1<<dim;
  }

  //! Return number of possible values for next position in empty multi index
  size_type subIndexCount() const
  {
    return gridView_.size(dim);
  }

  //! Return number possible values for next position in multi index
  size_type subIndexCount(const MultiIndex& index) const DUNE_FINAL
  {
    return gridView_.size(dim);
  }

  /**
   * \brief Return local view for basis
   *
   * Perhaps we must move the construction outside
   * of the global basis in order to calm the compiler
   * when instantiating the TMP constructing the local view.
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

protected:
  const GridView gridView_;
};


template<typename GV>
class PQ1NodalBasisLocalView
{
public:

  typedef PQ1NodalBasis<GV> GlobalBasis;
  typedef typename GlobalBasis::GridView GridView;
  typedef typename GlobalBasis::size_type size_type;
  typedef typename GlobalBasis::MultiIndex MultiIndex;
  typedef typename GridView::template Codim<0>::Entity Element;
  typedef PQ1NodalBasisLeafNode<GV> Tree;

  PQ1NodalBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis)
  {}

  //! Bind to element
  void bind(const Element& e)
  {
    element_ = &e;
    tree_.bind(e);
  }

  const Element& element() const
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  //! Unbind from element - only hint
  void unbind()
  {}

  /**
   * \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  const Element* element_;
  Tree tree_;
};


template<typename GV>
class PQ1NodalBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, GV::dimension, 1>::FiniteElementType,
    typename PQ1NodalBasis<GV>::size_type,
    typename PQ1NodalBasis<GV>::MultiIndex>
{
  typedef PQ1NodalBasis<GV> GlobalBasis;
  static const int dim = GV::dimension;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, 1> FiniteElementCache;
  typedef typename FiniteElementCache::FiniteElementType FE;
  typedef typename GlobalBasis::size_type ST;
  typedef typename GlobalBasis::MultiIndex MI;


  typedef typename GlobalBasis::LocalView LocalView;
  friend LocalView;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST,MI> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::MultiIndex MultiIndex;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  PQ1NodalBasisLeafNode(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    finiteElement_(0),
    element_(0)
  {}

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL // all nodes
  {
    return *element_;
  }

  const FiniteElement& finiteElement() const DUNE_FINAL // leaf nodes
  {
    return *finiteElement_;
  }

  //! size of subtree rooted in this node (element-local)
  size_type subTreeSize() const DUNE_FINAL // all nodes or leaf nodes only ?
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_->localBasis().size();
  }

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type maxSubTreeSize() const DUNE_FINAL // all nodes or leaf nodes only ?
  {
    return 1<<dim;
  }

  //! size of complete tree (element-local)
  size_type localSize() const DUNE_FINAL // all nodes
  {
    // We have localSize==subTreeSize because the tree consist of a single leaf node.
    return subTreeSize();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const DUNE_FINAL // all nodes
  {
    return i;
  }

  //! maximum size of complete tree for any element of the global basis
  size_type maxLocalSize() const DUNE_FINAL // all nodes
  {
    // We have maxLocalSize==maxSubTreeSize because the tree consist of a single leaf node.
    return maxSubTreeSize();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis (pair of multi-indices)
  const MultiIndex globalIndex(size_type i) const DUNE_FINAL // move to LocalView?
  {
    return { globalBasis_->gridView().indexSet().subIndex(
        *element_,
        finiteElement_->localCoefficients().localKey(i).subEntity(),
        dim) };
  }

  //! Generate multi indices for current subtree into range starting at it
  //! \param it iterator over a container of MultiIndex
  //! \return iterator past the last written element (STL-style)
  template<typename MultiIndexIterator>
  MultiIndexIterator generateMultiIndices(MultiIndexIterator it) const // move to LocalView?
  {
    size_type size = subTreeSize();
    for(size_type i=0; i<size; ++i)
    {
      (*it) = {globalBasis_->gridView().indexSet().subIndex(
        *element_,
        finiteElement_->localCoefficients().localKey(i).subEntity(),
        dim)};
      ++it;
    }
    return it;
  }

protected:

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
  }

  const GlobalBasis* globalBasis_;
  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH
