// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/std/final.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQ1NodeFactory
//   PQ1NodeIndexSet
//   PQ1Node
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, typename ST, typename TP>
class PQ1Node;

template<typename GV, class MI, class TP, class ST>
class PQ1NodeIndexSet;

template<typename GV, class MI, class ST>
class PQ1NodeFactory;

template<typename GV, class MI, class ST>
class PQ1NodeFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;

  template<class TP>
  using Node = PQ1Node<GV, size_type, TP>;

  template<class TP>
  using IndexSet = PQ1NodeIndexSet<GV, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

  /** \brief Constructor for a given grid view object */
  PQ1NodeFactory(const GridView& gv) :
    gridView_(gv)
  {}

  void initializeIndices()
  {}

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

  size_type size() const
  {
    return gridView_.size(dim);
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return size();
    if (prefix.size() == 1)
      return 0;
    DUNE_THROW(RangeError, "Method size() can only be called for prefixes of length up to one");
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return StaticPower<2,GV::dimension>::power;
  }

//protected:
  const GridView gridView_;
};



template<typename GV, typename ST, typename TP>
class PQ1Node :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, GV::dimension, 1>::FiniteElementType,
    ST,
    TP>
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<2,GV::dimension>::power;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, 1> FiniteElementCache;
  typedef typename FiniteElementCache::FiniteElementType FE;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST,TP> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;
  typedef typename Interface::TreePath TreePath;

  PQ1Node(const TreePath& treePath) :
    Interface(treePath),
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

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return size_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
    size_ = finiteElement_->size();
  }

protected:

  size_type size_;

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



template<typename GV, class MI, class TP, class ST>
class PQ1NodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQ1NodeFactory<GV, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  PQ1NodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    return {{ gridIndexSet.subIndex(element,localKey.subEntity(),dim) }};
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;
};

/** \brief Nodal basis of a scalar first-order Lagrangian finite element space
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam ST The type used for local indices; global indices are std::array<ST,1>
 */
template<typename GV, class ST = std::size_t>
using PQ1NodalBasis = DefaultGlobalBasis<PQ1NodeFactory<GV, std::array<ST, 1>, ST> >;

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1NODALBASIS_HH
