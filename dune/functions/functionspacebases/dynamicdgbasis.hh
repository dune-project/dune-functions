// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMICDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMICDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/dynamicordernode.hh>
#include <dune/functions/functionspacebases/qkdynamicordercache.hh>

#include <dune/typetree/treepath.hh>

namespace Dune {
namespace Functions {

/*
 * This is a dynamic order Qk node.
 * If one wants something different, one has to put a different type here.
 * This would be nicer, if it was a template argument of the PreBasis; but
 * one still needs the TP template...
 */
template<typename GV, typename TP>
using DynamicDGNode = DynamicOrderNode<GV, TP,
  DynamicQkFiniteElementCache<typename GV::ctype, double, GV::dimension>
>;

template<typename GV, class MI, class TP>
class DynamicDGNodeIndexSet;


/**
 * \brief PreBasis for a dynamic order DG basis.
 *
 * The actual finite elements are set through the DynamicDGNode type.
 * If one wants something different from Qk, one has to modify the node type.
 */
template<typename GV, class MI>
class DynamicDGPreBasis
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;
  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

  template<class TP>
  using Node = DynamicDGNode<GV, TP>;

  template<class TP>
  using IndexSet = DynamicDGNodeIndexSet<GV, MI, TP>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

  /** \brief Constructor for a given grid view object with constant
   * polynomial degree for all elements */
  DynamicDGPreBasis(const GridView& gv, int k=1) :
    gridView_(gv),
    mapper_(gv, Dune::mcmgElementLayout())
  {
    degrees_ = std::vector<int>(mapper_.size(), k);
    offSets_ = std::vector<size_type>(mapper_.size(),0);
  }

  /** \brief Constructor for a given grid view object
   * with a vector carrying the polymial degrees per element
   * (mapped according to the mapper scheme)
   */
  DynamicDGPreBasis(const GridView& gv, std::vector<int> degrees) :
    gridView_(gv),
    mapper_(gv, Dune::mcmgElementLayout()),
    degrees_(std::move(degrees))
  {
    offSets_ = std::vector<size_type>(mapper_.size(),0);
  }

  void initializeIndices()
  {
    // this is a first try: I don't know yet how reasonable this is
    // we use a node with a dummy treepath to get the sizes of the local finite elements
    // and compute a offset for each element. While we're at it, we also compute max size
    // and total size.
    auto node = Node<Dune::TypeTree::HybridTreePath<>>(Dune::TypeTree::HybridTreePath<>{}, &mapper_, &degrees_, &offSets_);
    size_type offSet = 0;

    maxNodeSize_=0;
    size_=0;

    for (const auto& e: elements(gridView_)) {
      node.bind(e);
      offSets_[mapper_.index(e)]=offSet;

      auto size = node.finiteElement().size();
      maxNodeSize_ = std::max(maxNodeSize_, static_cast<size_type>(size));

      offSet += size;
      size_ += size;
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  void update(const GridView& gv)
  {
    // This can not be used as updating the gridview without updating the degrees and will lead to undefined behaviour.
    DUNE_THROW(Dune::Exception, "update(gridView) cannot be used for a dynamic DG basis. Please recreate the whole basis object!");
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp, &mapper_, &degrees_, &offSets_};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    return size_;
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
    return maxNodeSize_;
  }

  /**
   * \brief Access to the mapper the basis uses to associate
   * polynomial degrees to elements.
   */
  const Mapper& mapper() const {
    return mapper_;
  }

  /**
   * \brief Get a mapper to associate
   * polynomial degrees to elements.
   */
  static Mapper mapper(const GV& gv) {
    return Mapper(gv, Dune::mcmgElementLayout());
  }

  /**
   * \brief Get a vector with the current element degrees mapped
   * according to mapper()
   */
  const auto& degrees() const {
    return degrees;
  }

//protected:
  GridView gridView_;

  Mapper mapper_;
  std::vector<int> degrees_;
  std::vector<size_type> offSets_;

  size_type size_;
  size_type maxNodeSize_;

}; // end class DynamicDGPreBasis



template<typename GV, class MI, class TP>
class DynamicDGNodeIndexSet
{
  // Cannot be an enum -- otherwise the switch statement below produces compiler warnings
  static const int dim = GV::dimension;

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = DynamicDGPreBasis<GV, MI>;

  using Node = typename PreBasis::template Node<TP>;

  DynamicDGNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis)
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
  template<typename It>
  It indices(It it) const
  {
    auto offSet = node_->offSet();

    for (size_type i = 0, end = size(); i < end ; ++i, ++it) {
      *it = MI{offSet + i};
    }
    return it;
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;
}; // end class DynamicDGNodeIndexSet


/** \brief Basis of a Qk DG finite element space with dynamic, possibly non-uniform local order.
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using DynamicQkDGBasis = DefaultGlobalBasis<DynamicDGPreBasis<GV, FlatMultiIndex<std::size_t>> >;



} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMICDGBASIS_HH
