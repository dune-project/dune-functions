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
//#include <dune/functions/functionspacebases/pq2nodalbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>

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

  using PQMultiIndex = std::array<size_type, 1>;
  using PQ1Factory = PQkNodeFactory<GV,1,PQMultiIndex>;
  using PQ2Factory = PQkNodeFactory<GV,2,PQMultiIndex>;

  using PQ1LocalIndexFunctor=ShiftedIdentity<size_type>;
  using PQ2LocalIndexFunctor=ShiftedIdentityWithStride<size_type>;

  using PQ1NodeIndexSet = typename PQ1Factory::template IndexSet<PQ1LocalIndexFunctor>;
  using PQ2NodeIndexSet = typename PQ2Factory::template IndexSet<PQ2LocalIndexFunctor>;

  TaylorHoodLocalIndexSet(const TaylorHoodIndexSet<GV> & indexSet) :
    pq1NodeIndexSet_(indexSet.basis_->pq1Factory_.template indexSet<PQ1LocalIndexFunctor>()),
    pq2NodeIndexSet_(indexSet.basis_->pq2Factory_.template indexSet<PQ2LocalIndexFunctor>())
  {}

  void bind(const TaylorHoodBasisLocalView<GV>& localView)
  {
    localView_ = & localView;

    pq1NodeIndexSet_.bind(localView_->tree().template child<1>());
    pq2NodeIndexSet_.bind(localView_->tree().template child<0>().child(0));
  }
  void unbind()
  {
    localView_ = nullptr;
    pq1NodeIndexSet_.unbind();
    pq2NodeIndexSet_.unbind();
  }

  size_type size() const
  {
    return dim*pq2NodeIndexSet_.size() + pq1NodeIndexSet_.size();
  }

  MultiIndex index(size_type localIndex) const
  {
    MultiIndex mi;
    size_type v_size = dim * pq2NodeIndexSet_.size();
    mi[0] = localIndex / v_size;
    if (mi[0] == 0)
    {
      size_type v_comp = localIndex % dim;
      size_type v_localIndex = localIndex / dim;
      mi[1] = pq2NodeIndexSet_.index(v_localIndex)[0] * dim + v_comp;
    }
    if (mi[0] == 1)
      mi[1] = pq1NodeIndexSet_.index(localIndex-v_size)[0];
    return mi;
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

private:
  const LocalView* localView_;
  PQ1NodeIndexSet pq1NodeIndexSet_;
  PQ2NodeIndexSet pq2NodeIndexSet_;
};

template<typename GV>
class TaylorHoodIndexSet
{
  static const int dim = GV::dimension;

  /** \brief The global FE basis that this is a view on */
  typedef TaylorHoodBasis<GV> GlobalBasis;

  TaylorHoodIndexSet(const GlobalBasis & basis) :
    basis_(&basis)
  {}

public:

  typedef TaylorHoodLocalIndexSet<GV> LocalIndexSet;

  typedef std::size_t size_type;

  /** \todo This enum has been added to the interface without prior discussion. */
  enum { multiIndexMaxSize = 2 };
  typedef std::array<size_type,2> MultiIndex;

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return dim * basis_->pq2Factory_.size()
      + basis_->pq1Factory_.size();
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
        return dim * basis_->pq2Factory_.size();
      if (prefix[0] == 1)
        return basis_->pq2Factory_.size();
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

  const GlobalBasis* basis_;
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
public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

protected:

  static const int dim = GV::dimension;

  using PQMultiIndex = std::array<size_type, 1>;
  using PQ1Factory = PQkNodeFactory<GV,1,PQMultiIndex>;
  using PQ2Factory = PQkNodeFactory<GV,2,PQMultiIndex>;

  using PQ1LocalIndexFunctor=ShiftedIdentity<size_type>;
  using PQ2LocalIndexFunctor=ShiftedIdentityWithStride<size_type>;

  using PQ2Node=PQkNodalBasisLeafNode<GV,2, PQ2LocalIndexFunctor >;

public:

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef TaylorHoodBasisLocalView<GV> LocalView;

  /** \brief Constructor for a given grid view object */
  TaylorHoodBasis(const GridView& gv) :
    pq1Factory_(gv),
    pq2Factory_(gv)
  {
    pq1Factory_.initializeIndices();
    pq2Factory_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return pq1Factory_.gridView();
  }

  TaylorHoodIndexSet<GV> indexSet() const
  {
    return TaylorHoodIndexSet<GV>(*this);
  }

  /**
   * \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

//private:
  friend TaylorHoodIndexSet<GV>;
  friend TaylorHoodBasisLocalView<GV>;

  PQ1Factory pq1Factory_;
  PQ2Factory pq2Factory_;
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
    tree_()
  {
  }

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = e;

    for(int i=0; i<GV::dimension; ++i)
      tree_.template child<0>().child(i).bind(element_);
    tree_.template child<1>().bind(element_);

//    size_type pressureOffset = tree_.template child<0>().size();
    size_type pressureOffset = tree_.template child<0>().child(0).size() * dim;

    tree_.template child<1>().localIndexFunctor() = ShiftedIdentity<std::size_t>(pressureOffset);

    // TODO: Implement this using generic tree-leaf loop
    size_ = dim*tree_.template child<0>().child(0).size() + tree_.template child<1>().size();
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    return element_;
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   * And indeed, in the TaylorHoodBasisView implementation this method does nothing.
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

  Tree& tree()
  {
    return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
    return size_;
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
    return dim*globalBasis_->pq2Factory_.maxNodeSize() + globalBasis_->pq1Factory_.maxNodeSize();
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
  Element element_;
  Tree tree_;
  size_type size_;
};

template<typename GV>
class TaylorHoodVelocityTree :
    public TypeTree::PowerNode<PQkNodalBasisLeafNode<GV,2, ShiftedIdentityWithStride<std::size_t> >, GV::dimension>
{
//  friend TaylorHoodBasisLocalView<GV>;

  using LocalIndexFunctor=ShiftedIdentityWithStride<std::size_t>;
  using Q2Node=PQkNodalBasisLeafNode<GV,2, LocalIndexFunctor >;
  using Base=TypeTree::PowerNode<Q2Node, GV::dimension>;

public:
  TaylorHoodVelocityTree() :
    Base()
  {
    for(int i=0; i<GV::dimension; ++i)
      this->setChild(i, std::make_shared<Q2Node>(LocalIndexFunctor(i, GV::dimension)));
  }
};

template<typename GV>
class TaylorHoodBasisTree :
    public TypeTree::CompositeNode<TaylorHoodVelocityTree<GV>,
                                   PQkNodalBasisLeafNode<GV,1,ShiftedIdentity<std::size_t>>>
{
  friend TaylorHoodBasisLocalView<GV>;

  using VelocityNode=TaylorHoodVelocityTree<GV>;
  using PressureNode=PQkNodalBasisLeafNode<GV,1,ShiftedIdentity<std::size_t>>;

  using Base=TypeTree::CompositeNode<VelocityNode, PressureNode>;

public:
  TaylorHoodBasisTree():
    Base()
  {
      this->template setChild<0>(std::make_shared<VelocityNode>());
      // 0 is just a dummy here. We must set the correct offset during bind (???).
      this->template setChild<1>(std::make_shared<PressureNode>(ShiftedIdentity<std::size_t>(0)));
  }
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
