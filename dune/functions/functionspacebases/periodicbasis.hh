// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PERIODICBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PERIODICBASIS_HH

#include <set>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune::Functions {

// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PeriodicPreBasis
//   PeriodicNodeIndexSet
//   PeriodicNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV>
class PeriodicNode;

template<typename GV, class MI>
class PeriodicNodeIndexSet;

template<typename GV, class MI>
class PeriodicPreBasis
{
  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  using Node = PeriodicNode<GV>;

  using IndexSet = PeriodicNodeIndexSet<GV, MI>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object
  PeriodicPreBasis(const GridView& gv) :
    gridView_(gv)
  {
    // Set up the index equivalence classes.
    // Without further identifications, each index is only equivalent to itself.
    indexEquivalenceClasses_.resize(gridView_.indexSet().size(dim));

    for (std::size_t i=0; i<indexEquivalenceClasses_.size(); i++)
      indexEquivalenceClasses_[i].insert(MI({i}));
  }

  void unifyIndexPair(MI a, MI b)
  {
    indexEquivalenceClasses_[a].insert(indexEquivalenceClasses_[b].begin(), indexEquivalenceClasses_[b].end());
    indexEquivalenceClasses_[b].insert(indexEquivalenceClasses_[a].begin(), indexEquivalenceClasses_[a].end());
  }

  void initializeIndices()
  {
    MI invalid = {std::numeric_limits<unsigned int>::max()};

    mappedIdx_.resize(indexEquivalenceClasses_.size());
    std::fill(mappedIdx_.begin(), mappedIdx_.end(), invalid);

    numIndices_ = 0;

    for (std::size_t i=0; i<indexEquivalenceClasses_.size(); i++)
    {
      auto someRepresentative = *indexEquivalenceClasses_[i].begin();

      if (mappedIdx_[someRepresentative] == invalid)
      {
        for (auto hostIdx : indexEquivalenceClasses_[i])
          mappedIdx_[hostIdx] = MI({numIndices_});
        numIndices_++;
      }
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{};
  }

  /**
   * \brief Create tree node index set
   *
   * Create an index set suitable for the tree node obtained
   * by makeNode().
   */
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }

  size_type size() const
  {
    return numIndices_;
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

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return numIndices_;
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return StaticPower<2,GV::dimension>::power;
  }

//protected:
  // For each idx of the host basis, mappedIdx_[idx] contains the corresponding
  // index in the index set with unified degrees of freedom
  std::vector<MI> mappedIdx_;

  std::size_t numIndices_;

  std::vector<std::set<MultiIndex> > indexEquivalenceClasses_;

  const GridView gridView_;
};



template<typename GV>
class PeriodicNode :
  public Functions::LeafBasisNode
{
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<2,GV::dimension>::power;

  using FiniteElementCache = typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, 1>;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  PeriodicNode() :
    finiteElement_(nullptr),
    element_(nullptr)
  {}

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
    finiteElement_ = &(cache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



template<typename GV, class MI>
class PeriodicNodeIndexSet
{
  enum {dim = GV::dimension};

public:

  using size_type = std::size_t;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using PreBasis = PeriodicPreBasis<GV, MI>;

  using Node = PeriodicNode<GV>;

  PeriodicNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis),
    node_(nullptr)
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
    assert(node_ != nullptr);
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = preBasis_->gridView().indexSet();
    const auto& element = node_->element();

    // Index on the host grid
    MultiIndex hostIdx{{gridIndexSet.subIndex(element,localKey.subEntity(),dim) }};

    return preBasis_->mappedIdx_[hostIdx];
  }

protected:
  const PreBasis* preBasis_;

  const Node* node_;
};



namespace BasisFactory {

namespace Impl {

class PeriodicPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  // \brief Default constructor
  PeriodicPreBasisFactory()
  {}

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return PeriodicPreBasis<GridView, MultiIndex>(gridView);
  }
};

} // end namespace BasisFactory::Impl



/**
 * \brief Create a pre-basis factory that can create a periodic pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
auto periodic()
{
  return Impl::PeriodicPreBasisFactory();
}

} // end namespace BasisFactory


/** \brief Meta basis that allows to identify degrees of freedom of a given other basis
 *
 * Main intended use are spaces with periodic boundary conditions.
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using PeriodicBasis = Functions::DefaultGlobalBasis<PeriodicPreBasis<GV, Functions::FlatMultiIndex<std::size_t> > >;

} // end namespace Dune::Functions

#endif // DUNE_FUFEM_PERIODICBASIS_HH
