// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/std/final.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {

template<typename GV>
class PQ1NodalBasisLocalView;

template<typename GV>
class PQ1NodalBasisLeafNode;

template<typename GV>
class PQ1IndexSet;

template<typename GV>
class PQ1LocalIndexSet
{
  static const int dim = GV::dimension;

public:
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef PQ1NodalBasisLocalView<GV> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  PQ1LocalIndexSet()
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const PQ1NodalBasisLocalView<GV>& localView)
  {
    localView_ = &localView;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return localView_->tree().finiteElement_->size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis (pair of multi-indices)
  const MultiIndex index(size_type i) const
  {
    return {{ localView_->globalBasis_->gridView().indexSet().subIndex(
        *(localView_->element_),
        localView_->tree().finiteElement_->localCoefficients().localKey(i).subEntity(),
        dim) }};
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

  const PQ1NodalBasisLocalView<GV>* localView_;
};

template<typename GV>
class PQ1IndexSet
{
  static const int dim = GV::dimension;

public:

  typedef PQ1LocalIndexSet<GV> LocalIndexSet;

  PQ1IndexSet(const GV& gridView)
    : size_(gridView.size(dim))
  {}

  std::size_t size() const
  {
    return size_;
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet();
  }

private:
  std::size_t size_;
};



/** \brief Nodal basis of a scalar first-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV>
class PQ1NodalBasis
: public GridViewFunctionSpaceBasis<GV,
                                    PQ1NodalBasisLocalView<GV>,
                                    PQ1IndexSet<GV>,
                                    std::array<std::size_t, 1> >
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
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

  PQ1IndexSet<GV> indexSet() const
  {
    return PQ1IndexSet<GV>(gridView_);
  }

  /**
   * \brief Return local view for basis
   *
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
  static const int dim = GV::dimension;

  // allow acces to the global function space
  friend PQ1LocalIndexSet<GV>;

public:

  typedef PQ1NodalBasis<GV> GlobalBasis;
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
   * In the case of a P1 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  typedef PQ1NodalBasisLeafNode<GV> Tree;

  /** \brief Construct local view for a given global finite element basis */
  PQ1NodalBasisLocalView(const GlobalBasis* globalBasis) :
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

  /**
   * \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /**
   * \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  Tree& tree()
  {
    return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return tree_.finiteElement_->size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * The method returns 3^dim, which is the number of degrees of freedom you get for cubes.
   */
  size_type maxSize() const
  {
    return 1<<dim;
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
class PQ1NodalBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, GV::dimension, 1>::FiniteElementType,
    typename PQ1NodalBasis<GV>::size_type>
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
  friend class PQ1LocalIndexSet<GV>;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  PQ1NodalBasisLeafNode(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    finiteElement_(nullptr),
    element_(nullptr)
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

  //! maximum size of subtree rooted in this node for any element of the global basis
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_->size();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const DUNE_FINAL
  {
    return i;
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
