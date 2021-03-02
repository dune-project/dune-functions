// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH


#include <tuple>
#include <optional>

#include <dune/common/concept.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/functions/functionspacebases/concepts.hh>



namespace Dune {
namespace Functions {

namespace Impl {

// Concept for a PreBasis implementing the indices() method.
//
// This concept is deprecated.
template<typename PreBasis>
using PreBasisHasIndicesConcept = decltype(std::declval<PreBasis>().indices(std::declval<typename PreBasis::Node>(), std::declval<std::vector<typename PreBasis::MultiIndex>>().begin()), true);

// Concept checking if a PreBasis implements the indices() method.
//
// This check is deprecated.
template<typename PreBasis>
using preBasisHasIndices = Std::is_detected<PreBasisHasIndicesConcept, PreBasis>;

// Backward compatible version of preBasis.indices(node, it)
// If the preBasis implements the method it's called directly.
// Otherwise this is forwarded to the old dprecated interface
// by creating a temporary NodeIndexSet. This may be expensive
// if the NodeIndexSet's member variables or bind() method are
// non-trivial.
//
// Once the old interface is gone and the new one is mandatory,
// this can be replaced by preBasis.indices(node, it).
template<typename PreBasis, typename Node, typename It>
It preBasisIndices(const PreBasis& preBasis, const Node& node, It multiIndices)
{
  if constexpr (preBasisHasIndices<PreBasis>{})
    return preBasis.indices(node, multiIndices);
  else
  {
    auto indexSet = preBasis().makeIndexSet();
    indexSet.bind(node);
    return indexSet.indices(multiIndices);
  }
}

// This class exists for backward compatibility only.
// A PreBasis providing indices(node,iterator) can
// simpe export a DefaultNodeIndexSet<PreBasis>
// which will forward the old-style interface to the
// new one.
//
// This class is deprecated.
template<class PB>
class DefaultNodeIndexSet
{
public:

  using size_type = std::size_t;
  using PreBasis = PB;
  using MultiIndex = typename PreBasis::MultiIndex;
  using Node = typename PreBasis::Node;

  DefaultNodeIndexSet(const PreBasis& preBasis) :
    preBasis_(&preBasis),
    node_(nullptr)
  {}

  void bind(const Node& node) {
    node_ = &node;
  }

  void unbind() {
    node_ = nullptr;
  }

  size_type size() const {
    assert(node_ != nullptr);
    return node_->size();
  }

  template<typename It>
  It indices(It it) const {
    assert(node_ != nullptr);
    return preBasis_->indices(*node_, it);
  }

protected:
  const PreBasis* preBasis_;
  const Node* node_;
};

}



/** \brief The restriction of a finite element basis to a single element */
template<class GB>
class DefaultLocalView
{
public:

  //! The global FE basis that this is a view on
  using GlobalBasis = GB;

  //! The grid view the global FE basis lives on
  using GridView = typename GlobalBasis::GridView;

  //! Type of the grid element we are bound to
  using Element = typename GridView::template Codim<0>::Entity;

  //! The type used for sizes
  using size_type = std::size_t;

  //! Tree of local finite elements / local shape function sets
  using Tree = typename GlobalBasis::PreBasis::Node;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = typename GlobalBasis::PreBasis::MultiIndex;

private:

  // The following helpers should be removed after 2.8
  template<typename NodeIndexSet_>
  using hasIndices = decltype(std::declval<NodeIndexSet_>().indices(std::declval<std::vector<typename NodeIndexSet_::MultiIndex>>().begin()));

  // A dummy NodeIndexSet to be used if the PreBasis provides
  // the new indices() method.
  struct DummyNodeIndexSet {};

  static auto makeIndexSet(const typename GlobalBasis::PreBasis& preBasis)
  {
    if constexpr (Impl::preBasisHasIndices<typename GlobalBasis::PreBasis>{})
      return DummyNodeIndexSet{};
    else
      return preBasis.makeIndexSet();
  }

  using NodeIndexSet = std::decay_t<decltype(makeIndexSet(std::declval<typename GlobalBasis::PreBasis>()))>;

public:

  /** \brief Construct local view for a given global finite element basis */
  DefaultLocalView(const GlobalBasis& globalBasis) :
    globalBasis_(&globalBasis),
    tree_(globalBasis_->preBasis().makeNode()),
    nodeIndexSet_(makeIndexSet(globalBasis_->preBasis()))
  {
    static_assert(models<Concept::BasisTree<GridView>, Tree>(), "Tree type passed to DefaultLocalView does not model the BasisNode concept.");
    initializeTree(tree_);
  }

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = e;
    bindTree(tree_, *element_);
    indices_.resize(size());
    if constexpr (Impl::preBasisHasIndices<typename GlobalBasis::PreBasis>{})
      globalBasis_->preBasis().indices(tree_, indices_.begin());
    else
    {
      nodeIndexSet_.bind(tree_);
      if constexpr (Std::is_detected<hasIndices,NodeIndexSet>{})
        nodeIndexSet_.indices(indices_.begin());
      else
        for (size_type i = 0 ; i < this->size() ; ++i)
          indices_[i] = nodeIndexSet_.index(i);
    }
  }

  /** \brief Return if the view is bound to a grid element
   */
  bool isBound() const {
    return static_cast<bool>(element_);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   */
  void unbind()
  {
    if constexpr (not Impl::preBasisHasIndices<typename GlobalBasis::PreBasis>{})
      nodeIndexSet_.unbind();
    element_.reset();
  }

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  /** \brief Total number of degrees of freedom on this element
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
   */
  size_type maxSize() const
  {
    return globalBasis_->preBasis().maxNodeSize();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    return indices_[i];
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

  const DefaultLocalView& rootLocalView() const
  {
    return *this;
  }

protected:
  const GlobalBasis* globalBasis_;
  std::optional<Element> element_;
  Tree tree_;
  NodeIndexSet nodeIndexSet_;
  std::vector<MultiIndex> indices_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DEFAULTLOCALVIEW_HH
