// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASCUBEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASCUBEBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/raviartthomas.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>

#include <dune/functions/pdelab-copy/raviartthomasfem.hh>

namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   RaviartThomasCubeNodeFactory
//   RaviartThomasCubeNodeIndexSet
//   RaviartThomasCubeNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename ST, typename TP>
class RaviartThomasCubeNode;

template<typename GV, int k, class MI, class TP, class ST>
class RaviartThomasCubeNodeIndexSet;

template<typename GV, int k, class MI, class ST>
class RaviartThomasCubeNodeFactory;

template<typename GV, int k, class MI, class ST>
class RaviartThomasCubeNodeFactory
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = ST;

  // Precompute the number of dofs per entity type depending on the entity's codimension
  std::vector<int> dofsPerCodim {dim*k, k+1}; // holds only for 2D!

  template<class TP>
  using Node = RaviartThomasCubeNode<GV, k, size_type, TP>;

  template<class TP>
  using IndexSet = RaviartThomasCubeNodeIndexSet<GV, k, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  /** \brief Constructor for a given grid view object */
  RaviartThomasCubeNodeFactory(const GridView& gv) :
    gridView_(gv)
  { }

  void initializeIndices()
  {
    CodimOffset_.resize(2);
    CodimOffset_[0] = 0;
    CodimOffset_[1] = CodimOffset_[0] + dofsPerCodim[0] * gridView_.size(0);
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
    return Node<TP>{tp, gridView_};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    return dofsPerCodim[0] * gridView_.size(0) + dofsPerCodim[1] * gridView_.size(1);
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return StaticPower<(k+1),GV::dimension>::power;
  }

//protected:
  const GridView gridView_;
  std::vector<size_t> CodimOffset_;

};



template<typename GV, int k, typename ST, typename TP>
class RaviartThomasCubeNode :
  public LeafBasisNode<ST, TP>
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<(k+1),GV::dimension>::power;

  using Base = LeafBasisNode<ST,TP>;

public:

  using size_type = ST;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElementMap = typename PDELab::detail::RaviartThomasLocalFiniteElementMapBaseSelector<GV, dim, Dune::GeometryType::cube, typename GV::ctype, double, k>::type;
  using FiniteElement =typename FiniteElementMap::Traits::FiniteElement;

  RaviartThomasCubeNode(const TreePath& treePath, const GV& gv) :
    Base(treePath),
    finiteElement_(nullptr),
    element_(nullptr),
    finiteElementMap_(gv)
  { }

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(finiteElementMap_.find(*element_));
    this->setSize(finiteElement_->size());
  }

protected:

  const FiniteElement* finiteElement_;
  const Element* element_;
  FiniteElementMap finiteElementMap_;
};



template<typename GV, int k, class MI, class TP, class ST>
class RaviartThomasCubeNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = RaviartThomasCubeNodeFactory<GV, k, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  RaviartThomasCubeNodeIndexSet(const NodeFactory& nodeFactory) :
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
  //! Assume dim \in \lbrace 2, 3 \rbrace.
  MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    size_t subentity = localKey.subEntity();
    size_t codim = localKey.codim();

    // throw if Element is no cube
    if (not(element.type().isCube())) DUNE_THROW(Dune::NotImplemented, "RaviartThomasCubeNodalBasis only implemented for cube elements.");

    if (codim==0) { // element dof
      return { nodeFactory_->CodimOffset_[codim] +
               nodeFactory_->dofsPerCodim[codim] * gridIndexSet.subIndex(element, subentity, codim) + localKey.index() };
    }

    // treat edges individually due to orientation
    if (codim==1) { // edge/face dof
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        // we have to reverse the numbering if the local triangle edge is
        // not aligned with the global edge
        size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(subentity,codim,0,dim),dim);
        size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(subentity,codim,1,dim),dim);
        bool flip = (v0 > v1);
        return { (flip)
          ? nodeFactory_->CodimOffset_[codim]
          + nodeFactory_->dofsPerCodim[codim] * gridIndexSet.subIndex(element,subentity,codim) + (nodeFactory_->dofsPerCodim[codim]-1) - localKey.index()
              : nodeFactory_->CodimOffset_[codim]
              + nodeFactory_->dofsPerCodim[codim] * gridIndexSet.subIndex(element,subentity,codim) + localKey.index() };
    }
    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the RaviartThomasCubeBasis");
  }

protected:
  const NodeFactory* nodeFactory_;
  const Node* node_;
};



namespace BasisBuilder {

namespace Imp {

template<std::size_t k>
struct RaviartThomasCubeNodeFactoryBuilder
{
  static const std::size_t requiredMultiIndexSize=1;

  template<class MultiIndex, class GridView, class size_type=std::size_t>
  auto build(const GridView& gridView)
    -> RaviartThomasCubeNodeFactory<GridView, k, MultiIndex, size_type>
  {
    return {gridView};
  }
};

} // end namespace BasisBuilder::Imp

template<std::size_t k>
Imp::RaviartThomasCubeNodeFactoryBuilder<k> rtcube()
{
  return{};
}

} // end namespace BasisBuilder



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Raviart Thomas finite element space on simplicial elements
 *
 * TODO
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k, class ST = std::size_t>
using RaviartThomasCubeNodalBasis = DefaultGlobalBasis<RaviartThomasCubeNodeFactory<GV, k, FlatMultiIndex<ST>, ST> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_RAVIARTTHOMASCUBEBASIS_HH
