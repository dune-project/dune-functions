// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/std/final.hh>

#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>

#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <dune/functions/functionspacebases/pq2nodalbasis.hh>

namespace Dune {
namespace Functions {


template<typename GV>
class TaylorHoodBasis;

template<typename GV>
class TaylorHoodBasisLocalView;

template<typename GV>
class TaylorHoodIndexSet;

template<typename GV>
class TaylorHoodVelocityTree;

template<typename GV>
class TaylorHoodBasisTree;

template<typename GV>
class TaylorHoodLocalIndexSet
{
  static const int dim = GV::dimension;

public:
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef TaylorHoodBasisLocalView<GV> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type,2> MultiIndex;

  TaylorHoodLocalIndexSet(const TaylorHoodIndexSet<GV> & indexSet)
    : indexSet_(indexSet)
    , pq1LocalIndexSet_(indexSet_.pq1IndexSet_.localIndexSet())
    , pq2LocalIndexSet_(indexSet_.pq2IndexSet_.localIndexSet())
  {}

  void bind(const TaylorHoodBasisLocalView<GV>& localView)
  {
    localView_ = & localView;
    pq1LocalIndexSet_.bind(localView.pq1localView_);
    pq2LocalIndexSet_.bind(localView.pq2localView_);
  }
  void unbind()
  {
    localView_ = nullptr;
    pq1LocalIndexSet_.unbind();
    pq2LocalIndexSet_.unbind();
  }

  size_type size() const
  {
    return dim*pq2LocalIndexSet_.size() + pq1LocalIndexSet_.size();
  }

  MultiIndex index(size_type localIndex) const
  {
    MultiIndex mi;
    size_type v_size = dim * pq2LocalIndexSet_.size();
    mi[0] = localIndex / v_size;
    if (mi[0] == 0)
    {
      size_type v_comp = (localIndex % v_size) / pq2LocalIndexSet_.size();
      size_type v_localIndex = (localIndex % v_size) % pq2LocalIndexSet_.size();
      mi[1] = v_comp * indexSet_.pq2IndexSet_.size() +
        pq2LocalIndexSet_.index(v_localIndex)[0];
    }
    if (mi[0] == 1)
    {
      size_type p_localIndex = localIndex % v_size;
      mi[1] = pq1LocalIndexSet_.index(p_localIndex)[0];
    }
    return mi;
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

private:
  const LocalView * localView_;
  TaylorHoodIndexSet<GV> indexSet_;
  typename PQ1NodalBasis<GV>::IndexSet::LocalIndexSet pq1LocalIndexSet_;
  typename PQ2NodalBasis<GV>::IndexSet::LocalIndexSet pq2LocalIndexSet_;
};

template<typename GV>
class TaylorHoodIndexSet
{
  static const int dim = GV::dimension;

  /** \brief The global FE basis that this is a view on */
  typedef TaylorHoodBasis<GV> GlobalBasis;

  TaylorHoodIndexSet(const GlobalBasis & basis)
    : pq1IndexSet_(basis.pq1nodalbasis_.indexSet())
    , pq2IndexSet_(basis.pq2nodalbasis_.indexSet())
  {}

public:

  typedef TaylorHoodLocalIndexSet<GV> LocalIndexSet;

  typedef std::size_t size_type;
#warning new interface enum
  enum { multiIndexMaxSize = 2 };
  typedef std::array<size_type,2> MultiIndex;

#warning new interface method
  size_type dimension() const
  {
    return dim * pq2IndexSet_.size()
      + pq1IndexSet_.size();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return 2;
  }

  //! Return number possible values for next position in multi index
  size_type size(Dune::ReservedVector<std::size_t, multiIndexMaxSize> prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix[0] == 0)
        return dim * pq2IndexSet_.size();
      if (prefix[0] == 1)
        return pq1IndexSet_.size();
    }
    assert(false);
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(*this);
  }

private:
  friend GlobalBasis;
  friend LocalIndexSet;

  typename PQ1NodalBasis<GV>::IndexSet pq1IndexSet_;
  typename PQ2NodalBasis<GV>::IndexSet pq2IndexSet_;
};



/** \brief Nodal basis of a scalar second-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV>
class TaylorHoodBasis
  : public GridViewFunctionSpaceBasis<GV,
                                      TaylorHoodBasisLocalView<GV>,
                                      TaylorHoodIndexSet<GV>,
                                      std::array<std::size_t, 2> >
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef TaylorHoodBasisLocalView<GV> LocalView;

  /** \brief Constructor for a given grid view object */
  TaylorHoodBasis(const GridView& gv)
    : pq1nodalbasis_(gv)
    , pq2nodalbasis_(gv)
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return pq1nodalbasis_.gridView();
  }

  TaylorHoodIndexSet<GV> indexSet() const
  {
    return TaylorHoodIndexSet<GV>(gridView());
  }

  /**
   * \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

private:
  friend TaylorHoodIndexSet<GV>;
  friend TaylorHoodBasisLocalView<GV>;

  PQ1NodalBasis<GV> pq1nodalbasis_;
  PQ2NodalBasis<GV> pq2nodalbasis_;
};


/** \brief The restriction of a finite element basis to a single element */
template<typename GV>
class TaylorHoodBasisLocalView
{
  static const int dim = GV::dimension;

public:
  /** \brief The global FE basis that this is a view on */
  typedef TaylorHoodBasis<GV> GlobalBasis;
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
  typedef TaylorHoodBasisTree<GV> Tree;

  /** \brief Construct local view for a given global finite element basis */
  TaylorHoodBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    pq1localView_(globalBasis->pq1nodalbasis_.localView()),
    pq2localView_(globalBasis->pq2nodalbasis_.localView()),
    velocityTree_(pq2localView_.tree()),
    tree_(velocityTree_, pq1localView_.tree())
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    pq1localView_.bind(e);
    pq2localView_.bind(e);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    pq1localView_.element();
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   * And indeed, in the TaylorHoodBasisView implementation this method does nothing.
   */
  void unbind()
  {
    pq1localView_.unbind();
    pq2localView_.unbind();
  }

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  Tree& tree()
  {
    return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
    return dim*pq2localView_.size() + pq1localView_.size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   */
  size_type maxSize() const
  {
    return dim*pq2localView_.maxSize() + pq1localView_.maxSize();
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  friend TaylorHoodLocalIndexSet<GV>;

  const GlobalBasis* globalBasis_;
  PQ1NodalBasisLocalView<GV> pq1localView_;
  PQ2NodalBasisLocalView<GV> pq2localView_;
  TaylorHoodVelocityTree<GV> velocityTree_;
  Tree tree_;
};

template<typename GV>
class TaylorHoodVelocityTree :
    public TypeTree::PowerNode<PQ2NodalBasisLeafNode<GV>, GV::dimension>
{
  friend TaylorHoodBasisLocalView<GV>;

  TaylorHoodVelocityTree(PQ2NodalBasisLeafNode<GV> & pq2leafNode) :
    TypeTree::PowerNode<PQ2NodalBasisLeafNode<GV>, GV::dimension>(
      pq2leafNode, false /* don't make copies */)
  {}
};

template<typename GV>
class TaylorHoodBasisTree :
    public TypeTree::CompositeNode<TaylorHoodVelocityTree<GV>,
                                   PQ1NodalBasisLeafNode<GV>>
{
  friend TaylorHoodBasisLocalView<GV>;

  typedef TypeTree::CompositeNode<TaylorHoodVelocityTree<GV>,
                                  PQ1NodalBasisLeafNode<GV>> Base;
  TaylorHoodBasisTree(TaylorHoodVelocityTree<GV> & velocityTree,
    PQ1NodalBasisLeafNode<GV> & pq1leafNode) :
    Base(velocityTree, pq1leafNode)
  {}
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
