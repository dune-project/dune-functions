// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ2NODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ2NODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {


template<typename GV>
class PQ2NodalBasisLocalView;

template<typename GV>
class PQ2NodalBasisLeafNode;



/** \brief Nodal basis of a scalar second-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV>
class PQ2NodalBasis
: public GridViewFunctionSpaceBasis<GV,
                                    PQ2NodalBasisLocalView<GV>,
                                    std::array<std::size_t, 1> >
{
  static const int dim = GV::dimension;

  template<int dim>
  struct P2MapperLayout
  {
    bool contains (Dune::GeometryType gt) const
    {
      // All hypercubes carry a degree of freedom (this includes vertices and edges)
      return gt.isCube();
    }
  };

  // Needs the mapper
  friend class PQ2NodalBasisLeafNode<GV>;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef PQ2NodalBasisLocalView<GV> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  /** \brief Constructor for a given grid view object */
  PQ2NodalBasis(const GridView& gv) :
    gridView_(gv),
    mapper_(gv)
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * max{GridViewLocalBasisView(e).tree().size() | e in GridView}
   *
   * The method returns 3^dim, which is the number of degrees of freedom you get for cubes.
   */
  size_type maxLocalSize() const DUNE_FINAL
  {
    return StaticPower<3,dim>::power;
  }

  //! Return the number of possible values for next position in empty multi index
  size_type subIndexCount() const
  {
    return mapper_.size();
  }

  //! Return number possible values for next position in multi index
  size_type subIndexCount(const MultiIndex& index) const DUNE_FINAL
  {
    return mapper_.size();
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

protected:
  const GridView gridView_;
  const MultipleCodimMultipleGeomTypeMapper<GridView, P2MapperLayout> mapper_;
};


/** \brief The restriction of a finite element basis to a single element */
template<typename GV>
class PQ2NodalBasisLocalView
{
public:
  /** \brief The global FE basis that this is a view on */
  typedef PQ2NodalBasis<GV> GlobalBasis;
  typedef typename GlobalBasis::GridView GridView;

  /** \brief The type used for sizes */
  typedef typename GlobalBasis::size_type size_type;

  /** \brief Type used to number the degrees of freedom
   *
   * In the case of mixed finite elements this really can be a multi-index, but for a standard
   * P2 space this is only a single-digit multi-index, i.e., it is an integer.
   */
  typedef typename GlobalBasis::MultiIndex MultiIndex;

  /** \brief Type of the grid element we are bound to */
  typedef typename GridView::template Codim<0>::Entity Element;

  /** \brief Tree of local finite elements / local shape function sets
   *
   * In the case of a P2 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  typedef PQ2NodalBasisLeafNode<GV> Tree;

  /** \brief Construct local view for a given global finite element basis */
  PQ2NodalBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = &e;
    tree_.bind(e);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   * And indeed, in the PQ2NodalBasisView implementation this method does nothing.
   */
  void unbind()
  {}

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /** \brief Return the global basis that we are a view on
   */
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
class PQ2NodalBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, GV::dimension, 2>::FiniteElementType,
    typename PQ2NodalBasis<GV>::size_type,
    typename PQ2NodalBasis<GV>::MultiIndex>
{
  typedef PQ2NodalBasis<GV> GlobalBasis;
  static const int dim = GV::dimension;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, 2> FiniteElementCache;
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

  PQ2NodalBasisLeafNode(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    finiteElement_(0),
    element_(0)
  {}

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
  {
    return *finiteElement_;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type subTreeSize() const DUNE_FINAL // all nodes or leaf nodes only ?
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    return finiteElement_->size();
#else
    return finiteElement_->localBasis().size();
#endif
  }

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type maxSubTreeSize() const DUNE_FINAL // all nodes or leaf nodes only ?
  {
    return StaticPower<3,dim>::power;
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
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    return { globalBasis_->mapper_.subIndex(
#else
    return { (size_t)globalBasis_->mapper_.map(
#endif
        *element_,
        finiteElement_->localCoefficients().localKey(i).subEntity(),
        finiteElement_->localCoefficients().localKey(i).codim()) };
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
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
      (*it) = {globalBasis_->mapper_.subIndex(
#else
      (*it) = {globalBasis_->mapper_.map(
#endif
        *element_,
        finiteElement_->localCoefficients().localKey(i).subEntity(),
        finiteElement_->localCoefficients().localKey(i).codim())};
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


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ2NODALBASIS_HH
