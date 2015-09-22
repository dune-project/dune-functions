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

template<typename GV, class ST, typename TP>
class TaylorHoodVelocityTree;

template<typename GV, class ST, typename TP>
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
  using Node = TaylorHoodBasisTree<GV, ST, TP>;

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
        return pq1Factory_.size();
    }
    if (prefix.size() == 1)
      return 0;
    assert(false);
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return dim * pq2Factory_.size() + pq1Factory_.size();
  }

  size_type maxNodeSize() const
  {
    return dim * pq2Factory_.maxNodeSize() + pq1Factory_.maxNodeSize();
  }

protected:
  const GridView gridView_;

  PQ1Factory pq1Factory_;
  PQ2Factory pq2Factory_;
};



template<typename GV, class ST, typename TP>
class TaylorHoodVelocityTree :
    public PowerBasisNode<ST, TP ,PQkNode<GV,2, ST, decltype(TypeTree::push_back(TP(), 0)) >, GV::dimension>
{
  using ComponentTreePath = decltype(TypeTree::push_back(TP(), 0));

  using PQ2Node = PQkNode<GV,2, ST, ComponentTreePath >;
  using Base = PowerBasisNode<ST, TP ,PQ2Node, GV::dimension>;

public:
  TaylorHoodVelocityTree(const TP& tp) :
    Base(tp)
  {
    for(int i=0; i<GV::dimension; ++i)
      this->setChild(i, std::make_shared<PQ2Node>(TypeTree::push_back(tp, i)));
  }
};

template<typename GV, class ST, typename TP>
class TaylorHoodBasisTree :
    public CompositeBasisNode<ST, TP,
      TaylorHoodVelocityTree<GV, ST, decltype(TypeTree::push_back<0>(TP()))>,
      PQkNode<GV,1,ST, decltype(TypeTree::push_back<1ul>(TP()))>
    >
{
  using VelocityTreePath = decltype(TypeTree::push_back<0ul>(TP()));
  using PressureTreePath = decltype(TypeTree::push_back<1ul>(TP()));

  using VelocityNode=TaylorHoodVelocityTree<GV, ST, VelocityTreePath>;
  using PressureNode=PQkNode<GV,1,ST, PressureTreePath>;

  using Base=CompositeBasisNode<ST, TP, VelocityNode, PressureNode>;

public:
  TaylorHoodBasisTree(const TP& tp):
    Base(tp)
  {
    using namespace Dune::TypeTree::Indices;
    this->template setChild<0>(std::make_shared<VelocityNode>(push_back(tp, _0)));
    this->template setChild<1>(std::make_shared<PressureNode>(push_back(tp, _1)));
  }
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



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar second-order Lagrangean finite element space
 *
 * \tparam GV The GridView that the space is defined on.
 */
template<typename GV, class ST = std::size_t>
using TaylorHoodBasis = DefaultGlobalBasis<TaylorHoodNodeFactory<GV, std::array<ST, 2>, ST> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
