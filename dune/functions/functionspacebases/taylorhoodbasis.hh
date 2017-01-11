// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>

#include <dune/functions/functionspacebases/nodes.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

namespace Dune {
namespace Functions {


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   TaylorHoodNodeFactory
//   TaylorHoodNodeIndexSet
//   TaylorHoodBasisTree
//   TaylorHoodVelocityTree
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, typename TP>
class TaylorHoodVelocityTree;

template<typename GV, typename TP>
class TaylorHoodBasisTree;

template<typename GV, class MI, class TP, bool HI>
class TaylorHoodNodeIndexSet;



/**
 * \brief Factory for lowest order Taylor-Hood basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The grid view that the FE basis is defined on
 * \tparam MI Type to be used for multi-indices
 * \tparam HI Flag to select hybrid indices
 *
 * \note This mainly serves as an example, since you can construct a factory with
 * the same functionality manually using
 * \code
 * static const int k = 1;
 * using VelocityFactory = PowerNodeFactory<MI,IMS,PQkNodeFactory<GV,k+1,MI>,dim>;
 * using PressureFactory = PQkNodeFactory<GV,k,MI>;
 * using TaylorHoodKNodeFactory = CompositeNodeFactory<MI, BlockedLexicographic, VelocityFactory, PressureFactory>;
 * \endcode
 * Where IMS is LeafBlockedInterleaved if HI is set and
 * FlatInterleaved otherwise.
 */
template<typename GV, class MI, bool HI=false>
class TaylorHoodNodeFactory
{
  static const bool useHybridIndices = HI;

  static const int dim = GV::dimension;

  template<class, class, class, bool>
  friend class TaylorHoodNodeIndexSet;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  template<class TP>
  using Node = TaylorHoodBasisTree<GV, TP>;

  //! Template mapping root tree path to type of created tree node index set
  template<class TP>
  using IndexSet = TaylorHoodNodeIndexSet<GV, MI, TP, HI>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 3>;

private:

  using PQMultiIndex = std::array<size_type, 1>;
  using PQ1Factory = PQkNodeFactory<GV,1,PQMultiIndex>;
  using PQ2Factory = PQkNodeFactory<GV,2,PQMultiIndex>;

public:

  //! Constructor for a given grid view object
  TaylorHoodNodeFactory(const GridView& gv) :
    gridView_(gv),
    pq1Factory_(gv),
    pq2Factory_(gv)
  {}

  //! Initialize the global indices
  void initializeIndices()
  {
    pq1Factory_.initializeIndices();
    pq2Factory_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    pq1Factory_.update(gv);
    pq2Factory_.update(gv);
  }

  /**
   * \brief Create tree node with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * By passing a non-trivial root tree path this can be used
   * to create a node suitable for being placed in a tree at
   * the position specified by the root tree path.
   */
  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp};
  }

  /**
   * \brief Create tree node index set with given root tree path
   *
   * \tparam TP Type of root tree path
   * \param tp Root tree path
   *
   * Create an index set suitable for the tree node obtained
   * by node(tp).
   */
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

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    return sizeImp<useHybridIndices>(prefix);
  }

private:

  template<bool hi,
    typename std::enable_if<not hi,int>::type = 0>
  size_type sizeImp(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix[0] == 0)
        return dim * pq2Factory_.size();
      if (prefix[0] == 1)
        return pq1Factory_.size();
    }
    if (prefix.size() == 2)
      return 0;
    assert(false);
  }

  template<bool hi,
    typename std::enable_if<hi,int>::type = 0>
  size_type sizeImp(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix[0] == 0)
        return pq2Factory_.size();
      if (prefix[0] == 1)
        return pq1Factory_.size();
    }
    if (prefix.size() == 2)
    {
      if (prefix[0] == 0)
        return dim;
      if (prefix[0] == 1)
        return 0;
    }
    if (prefix.size() == 3)
      return 0;
    assert(false);
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return dim * pq2Factory_.size() + pq1Factory_.size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return dim * pq2Factory_.maxNodeSize() + pq1Factory_.maxNodeSize();
  }

protected:
  GridView gridView_;

  PQ1Factory pq1Factory_;
  PQ2Factory pq2Factory_;
};



template<typename GV, typename TP>
class TaylorHoodVelocityTree :
    public PowerBasisNode<std::size_t, TP ,PQkNode<GV,2, decltype(TypeTree::push_back(TP(), 0)) >, GV::dimension>
{
  using ComponentTreePath = decltype(TypeTree::push_back(TP(), 0));

  using PQ2Node = PQkNode<GV,2, ComponentTreePath >;
  using Base = PowerBasisNode<std::size_t, TP ,PQ2Node, GV::dimension>;

public:
  TaylorHoodVelocityTree(const TP& tp) :
    Base(tp)
  {
    for(int i=0; i<GV::dimension; ++i)
      this->setChild(i, std::make_shared<PQ2Node>(TypeTree::push_back(tp, i)));
  }
};

template<typename GV, typename TP>
class TaylorHoodBasisTree :
    public CompositeBasisNode<std::size_t, TP,
      TaylorHoodVelocityTree<GV, decltype(TypeTree::push_back<0>(TP()))>,
      PQkNode<GV,1, decltype(TypeTree::push_back<1ul>(TP()))>
    >
{
  using VelocityTreePath = decltype(TypeTree::push_back<0ul>(TP()));
  using PressureTreePath = decltype(TypeTree::push_back<1ul>(TP()));

  using VelocityNode=TaylorHoodVelocityTree<GV, VelocityTreePath>;
  using PressureNode=PQkNode<GV,1, PressureTreePath>;

  using Base=CompositeBasisNode<std::size_t, TP, VelocityNode, PressureNode>;

public:
  TaylorHoodBasisTree(const TP& tp):
    Base(tp)
  {
    using namespace Dune::TypeTree::Indices;
    this->template setChild<0>(std::make_shared<VelocityNode>(push_back(tp, _0)));
    this->template setChild<1>(std::make_shared<PressureNode>(push_back(tp, _1)));
  }
};



template<typename GV, class MI, class TP, bool HI>
class TaylorHoodNodeIndexSet
{
  static const bool useHybridIndices = HI;

  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = TaylorHoodNodeFactory<GV, MI, HI>;

  using Node = typename NodeFactory::template Node<TP>;

  using PQ1TreePath = typename TypeTree::Child<Node,1>::TreePath;
  using PQ2TreePath = typename TypeTree::Child<Node,0,0>::TreePath;

  using PQ1NodeIndexSet = typename NodeFactory::PQ1Factory::template IndexSet<PQ1TreePath>;
  using PQ2NodeIndexSet = typename NodeFactory::PQ2Factory::template IndexSet<PQ2TreePath>;

  TaylorHoodNodeIndexSet(const NodeFactory & nodeFactory) :
    nodeFactory_(&nodeFactory),
    pq1NodeIndexSet_(nodeFactory_->pq1Factory_.template indexSet<PQ1TreePath>()),
    pq2NodeIndexSet_(nodeFactory_->pq2Factory_.template indexSet<PQ2TreePath>())
  {}

  void bind(const Node& node)
  {
    using namespace TypeTree::Indices;
    node_ = &node;
    pq1NodeIndexSet_.bind(node.child(_1));
    pq2NodeIndexSet_.bind(node.child(_0, 0));
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
  }

  MultiIndex index(size_type localIndex) const
  {
    return indexImp<useHybridIndices>(localIndex);
  }

  template<bool hi,
    typename std::enable_if<not hi,int>::type = 0>
  MultiIndex indexImp(size_type localIndex) const
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

  template<bool hi,
    typename std::enable_if<hi,int>::type = 0>
  MultiIndex indexImp(size_type localIndex) const
  {
    MultiIndex mi;

    size_type velocityComponentSize = pq2NodeIndexSet_.size();
    size_type pressureOffset = velocityComponentSize * dim;

    mi.push_back(localIndex / pressureOffset);
    if (mi[0] == 0)
    {
      size_type v_comp = localIndex / velocityComponentSize;
      size_type v_localIndex = localIndex % velocityComponentSize;
      mi.push_back(pq2NodeIndexSet_.index(v_localIndex)[0]);
      mi.push_back(v_comp);
    }
    if (mi[0] == 1)
      mi.push_back(pq1NodeIndexSet_.index(localIndex-pressureOffset)[0]);
    return mi;
  }

private:
  const NodeFactory* nodeFactory_;
  PQ1NodeIndexSet pq1NodeIndexSet_;
  PQ2NodeIndexSet pq2NodeIndexSet_;

  const Node* node_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/**
 * \brief Nodal basis for a lowest order Taylor-Hood Lagrangean finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on.
 *
 * \note This mainly serves as an example, since you can construct a basis with
 * the same functionality manually using
 * \code
 * static const int k = 1;
 * auto taylorHoodBasis = makeBasis(
 *   gridView,
 *   composite(
 *     power<dim>(
 *       lagrange<k+1>(),
 *       flatInterleaved()),
 *     lagrange<k>()
 *   ));
 * \endcode
 */
template<typename GV>
using TaylorHoodBasis = DefaultGlobalBasis<TaylorHoodNodeFactory<GV, std::array<std::size_t, 2>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
