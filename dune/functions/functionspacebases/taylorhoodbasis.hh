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

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/defaultglobalindexset.hh>

namespace Dune {
namespace Functions {


template<typename GV, class ST = std::size_t>
class TaylorHoodBasis;

template<typename GV>
class TaylorHoodBasisLocalView;

template<typename GV, typename TP>
class TaylorHoodVelocityTree;

template<typename GV, typename TP>
class TaylorHoodBasisTree;

template<typename GV, class MI, class TP, class ST>
class TaylorHoodNodeIndexSet;



template<typename GV, class MI, class ST>
class TaylorHoodNodeFactory
{
  static const int dim = GV::dimension;

  template<class, class, class, class>
  friend class TaylorHoodNodeIndexSet;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;


  template<class TP>
  using Node = TaylorHoodBasisTree<GV, TP>;

  template<class TP>
  using IndexSet = TaylorHoodNodeIndexSet<GV, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

private:

  using PQMultiIndex = std::array<size_type, 1>;
  using PQ1Factory = PQkNodeFactory<GV,1,PQMultiIndex,ST>;
  using PQ2Factory = PQkNodeFactory<GV,2,PQMultiIndex,ST>;

public:

  /** \brief Constructor for a given grid view object */
  TaylorHoodNodeFactory(const GridView& gv) :
    gridView_(gv),
    pq1Factory_(gv),
    pq2Factory_(gv)
  {}


  void initializeIndices()
  {
    pq1Factory_.initializeIndices();
    pq2Factory_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return 2;
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix[0] == 0)
        return dim * pq2Factory_.size();
      if (prefix[0] == 1)
        return pq2Factory_.size();
    }
    assert(false);
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return dim * pq2Factory_.size() + pq1Factory_.size();
  }

  size_type maxNodeSize() const
  {
  }

//protected:
  const GridView gridView_;

  PQ1Factory pq1Factory_;
  PQ2Factory pq2Factory_;
};



template<typename GV, class MI, class TP, class ST>
class TaylorHoodNodeIndexSet
{
  static const int dim = GV::dimension;

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = TaylorHoodNodeFactory<GV, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  using PQ1TreePath = typename ChildType<Node,1>::TreePath;
  using PQ2TreePath = typename ChildType<Node,0,0>::TreePath;

  using PQ1NodeIndexSet = typename NodeFactory::PQ1Factory::template IndexSet<PQ1TreePath>;
  using PQ2NodeIndexSet = typename NodeFactory::PQ2Factory::template IndexSet<PQ2TreePath>;

  TaylorHoodNodeIndexSet(const NodeFactory & nodeFactory) :
    nodeFactory_(&nodeFactory),
    pq1NodeIndexSet_(nodeFactory_->pq1Factory_.template indexSet<PQ1TreePath>()),
    pq2NodeIndexSet_(nodeFactory_->pq2Factory_.template indexSet<PQ2TreePath>())
  {}

  void bind(const Node& node)
  {
    node_ = &node;
    pq1NodeIndexSet_.bind(node_->template child<1>());
    pq2NodeIndexSet_.bind(node_->template child<0>().child(0));
  }

  void unbind()
  {
    node_ = nullptr;
    pq1NodeIndexSet_.unbind();
    pq2NodeIndexSet_.unbind();
  }

  size_type size() const
  {
    return node_->size();
//    return dim*pq2NodeIndexSet_.size() + pq1NodeIndexSet_.size();
  }

  MultiIndex index(size_type localIndex) const
  {
    MultiIndex mi;

    size_type velocityComponentSize = pq2NodeIndexSet_.size();
    size_type pressureOffset = velocityComponentSize * dim;

    mi[0] = localIndex / pressureOffset;
    if (mi[0] == 0)
    {
      size_type v_comp = localIndex / velocityComponentSize;
      size_type v_localIndex = localIndex % velocityComponentSize;
      mi[1] = pq2NodeIndexSet_.index(v_localIndex)[0] * dim + v_comp;
    }
    if (mi[0] == 1)
      mi[1] = pq1NodeIndexSet_.index(localIndex-pressureOffset)[0];
    return mi;
  }

private:
  const NodeFactory* nodeFactory_;
  PQ1NodeIndexSet pq1NodeIndexSet_;
  PQ2NodeIndexSet pq2NodeIndexSet_;

  const Node* node_;
};




/** \brief Nodal basis of a scalar second-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV, class ST>
class TaylorHoodBasis
  : public GridViewFunctionSpaceBasis<GV,
                                      TaylorHoodBasisLocalView<GV>,
                                      DefaultGlobalIndexSet<TaylorHoodBasisLocalView<GV>, TaylorHoodNodeFactory<GV, std::array<ST, 2>, ST> >,
                                      std::array<std::size_t, 2> >
{
public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;

  using size_type = ST;

protected:

  static const int dim = GV::dimension;

  using PQMultiIndex = std::array<size_type, 1>;
  using PQ1Factory = PQkNodeFactory<GV,1,PQMultiIndex,ST>;
  using PQ2Factory = PQkNodeFactory<GV,2,PQMultiIndex,ST>;

public:

  using MultiIndex = std::array<size_type, 2>;
  using NodeFactory = TaylorHoodNodeFactory<GV, MultiIndex, size_type>;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  using LocalView = TaylorHoodBasisLocalView<GV>;

  using GlobalIndexSet = DefaultGlobalIndexSet<LocalView, NodeFactory>;


  /** \brief Constructor for a given grid view object */
  TaylorHoodBasis(const GridView& gv) :
    nodeFactory_(gv)
  {
    nodeFactory_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return nodeFactory_.gridView();
  }

  GlobalIndexSet indexSet() const
  {
    return GlobalIndexSet(nodeFactory_);
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
//protected:
  friend TaylorHoodBasisLocalView<GV>;

  NodeFactory nodeFactory_;
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

  using TreePath = std::tuple<>;

  /** \brief Tree of local finite elements / local shape function sets
   *
   * In the case of a P2 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  typedef TaylorHoodBasisTree<GV,TreePath> Tree;

  /** \brief Construct local view for a given global finite element basis */
  TaylorHoodBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(TreePath())
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
    bindTree(tree_, element_, 0);
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
    return tree_.size();
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
    return dim*globalBasis_->p.nodeFactory_.q2Factory_.maxNodeSize() + globalBasis_->nodeFactory_.pq1Factory_.maxNodeSize();
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:

  const GlobalBasis* globalBasis_;
  Element element_;
  Tree tree_;
};

template<typename GV, typename TP>
class TaylorHoodVelocityTree :
    public PowerBasisNode<std::size_t, TP ,PQkNode<GV,2, std::size_t, decltype(extendTreePath(TP(), int())) >, GV::dimension>
{
  using ComponentTreePath = decltype(extendTreePath(TP(), int()));

  using PQ2Node = PQkNode<GV,2, std::size_t, ComponentTreePath >;
  using Base = PowerBasisNode<std::size_t, TP ,PQ2Node, GV::dimension>;

public:
  TaylorHoodVelocityTree(const TP& tp) :
    Base(tp)
  {
    for(int i=0; i<GV::dimension; ++i)
      this->setChild(i, std::make_shared<PQ2Node>(extendTreePath(tp, i)));
  }
};

template<typename GV, typename TP>
class TaylorHoodBasisTree :
    public CompositeBasisNode<std::size_t, TP,
      TaylorHoodVelocityTree<GV, decltype(extendTreePath<0>(TP()))>,
      PQkNode<GV,1,std::size_t, decltype(extendTreePath<1>(TP()))>
    >
{
  using VelocityTreePath = decltype(extendTreePath<0>(TP()));
  using PressureTreePath = decltype(extendTreePath<1>(TP()));


  using VelocityNode=TaylorHoodVelocityTree<GV, VelocityTreePath>;
  using PressureNode=PQkNode<GV,1,std::size_t, PressureTreePath>;

  using Base=CompositeBasisNode<std::size_t, TP, VelocityNode, PressureNode>;

public:
  TaylorHoodBasisTree(const TP& tp):
    Base(tp)
  {
    using namespace StaticIndices;

    this->template setChild<0>(std::make_shared<VelocityNode>(extendTreePath(tp, _0)));
    this->template setChild<1>(std::make_shared<PressureNode>(extendTreePath(tp, _1)));
  }
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
